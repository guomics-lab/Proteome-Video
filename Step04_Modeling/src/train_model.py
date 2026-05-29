# train_model.py — GAM_LatentODE training script (reference, pre-trained weights provided)

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import argparse
import functools
import json
from datetime import datetime
from tqdm import tqdm
from torchdiffeq import odeint_adjoint


# ── Model Architecture ───────────────────────────────────────────────────────
class EncoderRNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, latent_dim):
        super().__init__()
        self.gru = nn.GRU(input_dim, hidden_dim, batch_first=True)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x):
        _, hidden = self.gru(x)
        hidden = hidden.squeeze(0)
        return self.fc_mu(hidden), self.fc_logvar(hidden)


class Decoder(nn.Module):
    def __init__(self, latent_dim, output_dim, hidden_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
            nn.Linear(hidden_dim, output_dim))

    def forward(self, z):
        return self.net(z)


class ConditionalODEFunc(nn.Module):
    def __init__(self, latent_dim, hidden_dim, condition_dim):
        super().__init__()
        self.baseline_net = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim), nn.Tanh(),
            nn.Linear(hidden_dim, latent_dim))
        self.diff_net = nn.Sequential(
            nn.Linear(latent_dim + condition_dim, hidden_dim), nn.Tanh(),
            nn.Linear(hidden_dim, latent_dim))

    def forward(self, t, z, condition_embedding):
        cond_emb_expanded = condition_embedding.expand(z.shape[0], -1)
        return self.baseline_net(z) + self.diff_net(torch.cat([z, cond_emb_expanded], dim=1))


class GAM_LatentODE(nn.Module):
    def __init__(self, input_dim, latent_dim=64, hidden_dim=128, num_conditions=3, condition_dim=16):
        super().__init__()
        self.encoder = EncoderRNN(input_dim, hidden_dim, latent_dim)
        self.ode_func = ConditionalODEFunc(latent_dim, hidden_dim, condition_dim)
        self.decoder = Decoder(latent_dim, input_dim, hidden_dim)
        self.dynamic_embedding = nn.Embedding(num_conditions, condition_dim)
        self.shift_embedding = nn.Embedding(num_conditions, input_dim)
        self.dynamic_embedding.weight.data[0].zero_()
        self.shift_embedding.weight.data[0].zero_()
        for emb in [self.dynamic_embedding, self.shift_embedding]:
            emb.weight.register_hook(lambda g, idx=0: g * torch.tensor(
                [0 if i == idx else 1 for i in range(g.shape[0])], device=g.device).unsqueeze(1))

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x_obs, condition_id, time_to_predict):
        mu, logvar = self.encoder(x_obs)
        z0 = self.reparameterize(mu, logvar)
        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        dynamic_emb = self.dynamic_embedding(condition_id)
        ode_func = functools.partial(self.ode_func, condition_embedding=dynamic_emb)
        z_t = odeint_adjoint(ode_func, z0, time_to_predict.to(x_obs.device),
                             method='dopri5', adjoint_params=tuple(self.ode_func.parameters()))
        z_traj = z_t.permute(1, 0, 2)
        y_pred = self.decoder(z_traj) + self.shift_embedding(condition_id).unsqueeze(1)
        return y_pred, kl_loss


# ── Data Loading ─────────────────────────────────────────────────────────────
def load_raw_data(data_dir):
    C = torch.load(os.path.join(data_dir, "Fourier_timeseries_C.pt"), weights_only=False).T
    A = torch.load(os.path.join(data_dir, "Fourier_timeseries_A.pt"), weights_only=False).T
    H = torch.load(os.path.join(data_dir, "Fourier_timeseries_H.pt"), weights_only=False).T
    min_len = min(C.shape[0], A.shape[0], H.shape[0])
    full_data = torch.stack([C[:min_len], A[:min_len], H[:min_len]], dim=0)
    full_time = torch.arange(min_len, dtype=torch.float32)
    return full_data, full_time


# ── Training ─────────────────────────────────────────────────────────────────
def train(data_dir, output_dir, obs_len=12, latent_dim=128, hidden_dim=256,
          condition_dim=32, lr=5e-4, weight_decay=1e-5, epochs=5000,
          kl_anneal_epochs=1000, kl_weight_max=0.01):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    os.makedirs(output_dir, exist_ok=True)

    full_data, full_time = load_raw_data(data_dir)
    pred_len = full_data.shape[1] - obs_len
    observed_data = full_data[:, :obs_len, :]
    truth_future = full_data[:, obs_len:, :]
    time_future = full_time[obs_len:]

    input_dim = full_data.shape[2]
    num_conditions = full_data.shape[0]

    print(f"Data: {full_data.shape}  obs={obs_len}  device={device}")

    model = GAM_LatentODE(input_dim, latent_dim, hidden_dim, num_conditions, condition_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    mse = nn.MSELoss()

    best_loss = float('inf')
    history = {'epoch': [], 'total_loss': [], 'recon_loss': [], 'kl_loss': []}

    pbar = tqdm(range(epochs), desc="Training")
    for epoch in pbar:
        model.train()
        optimizer.zero_grad()

        kl_w = min(kl_weight_max, kl_weight_max * epoch / kl_anneal_epochs)
        total_loss_val, recon_val, kl_val = 0, 0, 0

        for i in range(num_conditions):
            x_obs = observed_data[i].unsqueeze(0).to(device)
            y_true = truth_future[i].unsqueeze(0).to(device)
            cond_id = torch.tensor([i], device=device)
            y_pred, kl_loss_batch = model(x_obs, cond_id, time_future)
            recon = mse(y_pred, y_true)
            loss = recon + kl_w * kl_loss_batch
            loss.backward()
            total_loss_val += loss.item()
            recon_val += recon.item()
            kl_val += kl_loss_batch.item()

        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        avg_loss = total_loss_val / num_conditions
        avg_recon = recon_val / num_conditions
        avg_kl = kl_val / num_conditions

        history['epoch'].append(epoch)
        history['total_loss'].append(avg_loss)
        history['recon_loss'].append(avg_recon)
        history['kl_loss'].append(avg_kl)

        pbar.set_postfix({"Loss": f"{avg_loss:.4f}", "Recon": f"{avg_recon:.4f}"})

        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'best_loss': best_loss,
                'model_params': {
                    'input_dim': input_dim, 'latent_dim': latent_dim,
                    'hidden_dim': hidden_dim, 'num_conditions': num_conditions,
                    'condition_dim': condition_dim},
                'timestamp': datetime.now().isoformat()},
                os.path.join(output_dir, "best_model.pth"))

    print(f"Best loss: {best_loss:.6f}")
    with open(os.path.join(output_dir, "training_history.json"), 'w') as f:
        json.dump(history, f, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default='./checkpoints')
    parser.add_argument('--obs_len', type=int, default=12)
    parser.add_argument('--epochs', type=int, default=5000)
    args = parser.parse_args()

    torch.manual_seed(42)
    np.random.seed(42)
    train(args.data_dir, args.output_dir, obs_len=args.obs_len, epochs=args.epochs)
