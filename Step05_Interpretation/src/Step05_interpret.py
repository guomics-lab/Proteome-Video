# Step05_interpret.py — Latent ODE inference, PCA, gradients; calls Step05_GSEA.R for figures

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os, sys, re, subprocess, functools, json, warnings
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
from torchdiffeq import odeint_adjoint

warnings.filterwarnings("ignore")

BASE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(BASE, "..", "output")
os.makedirs(OUTPUT, exist_ok=True)

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 12

SANDBOX = os.path.normpath(os.path.join(BASE, "..", ".."))

# ── Model Architecture ──────────────────────────────────────────────────────
class EncoderRNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, latent_dim):
        super().__init__()
        self.gru = nn.GRU(input_dim, hidden_dim, batch_first=True)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)
    def forward(self, x):
        _, h = self.gru(x); h = h.squeeze(0)
        return self.fc_mu(h), self.fc_logvar(h)

class Decoder(nn.Module):
    def __init__(self, latent_dim, output_dim, hidden_dim):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(latent_dim, hidden_dim), nn.ReLU(),
                                 nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
                                 nn.Linear(hidden_dim, output_dim))
    def forward(self, z): return self.net(z)

class ConditionalODEFunc(nn.Module):
    def __init__(self, latent_dim, hidden_dim, condition_dim):
        super().__init__()
        self.baseline_net = nn.Sequential(nn.Linear(latent_dim, hidden_dim), nn.Tanh(),
                                          nn.Linear(hidden_dim, latent_dim))
        self.diff_net = nn.Sequential(nn.Linear(latent_dim + condition_dim, hidden_dim), nn.Tanh(),
                                      nn.Linear(hidden_dim, latent_dim))
    def forward(self, t, z, cond_emb):
        return self.baseline_net(z) + self.diff_net(torch.cat([z, cond_emb.expand(z.shape[0], -1)], dim=1))

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
    def reparameterize(self, mu, logvar):
        return mu + torch.exp(0.5 * logvar) * torch.randn_like(logvar)
    def forward(self, x_obs, condition_id, time_to_predict):
        mu, logvar = self.encoder(x_obs)
        z0 = self.reparameterize(mu, logvar)
        kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        ode_f = functools.partial(self.ode_func, cond_emb=self.dynamic_embedding(condition_id))
        z_t = odeint_adjoint(ode_f, z0, time_to_predict.to(x_obs.device),
                             method='dopri5', adjoint_params=tuple(self.ode_func.parameters()))
        y = self.decoder(z_t.permute(1, 0, 2))
        return y + self.shift_embedding(condition_id).unsqueeze(1), kl
    def get_latent(self, x_obs, condition_id, time_points):
        mu, _ = self.encoder(x_obs)
        ode_f = functools.partial(self.ode_func, cond_emb=self.dynamic_embedding(condition_id))
        z_t = odeint_adjoint(ode_f, mu, time_points.to(x_obs.device),
                             method='dopri5', adjoint_params=tuple(self.ode_func.parameters()))
        return z_t.permute(1, 0, 2)


def preprocess_data():
    tsv_path = os.path.join(SANDBOX, "Step02_Clustering", "output", "downsampled_median.tsv")
    df = pd.read_csv(tsv_path, sep='\t', index_col=0)
    print(f"Loaded: {df.shape[0]} proteins x {df.shape[1]} samples")

    groups = {'C': [], 'A': [], 'H': []}
    for col in df.columns:
        m = re.match(r'([CAH])(\d+)', col)
        if m: groups[m.group(1)].append((col, int(m.group(2))))
    for k in groups: groups[k].sort(key=lambda x: x[1])

    for cond, cols in [('C', groups['C']), ('A', groups['A']), ('H', groups['H'])]:
        t = torch.tensor(df[[c[0] for c in cols]].values, dtype=torch.float32)
        torch.save(t, os.path.join(OUTPUT, f"Fourier_timeseries_{cond}.pt"))

    time_idx = [t for _, t in groups['C']]
    torch.save(torch.tensor(time_idx, dtype=torch.float32),
               os.path.join(OUTPUT, "Fourier_time_points.pt"))


def load_tensors():
    C = torch.load(os.path.join(OUTPUT, "Fourier_timeseries_C.pt"), weights_only=False).T
    A = torch.load(os.path.join(OUTPUT, "Fourier_timeseries_A.pt"), weights_only=False).T
    H = torch.load(os.path.join(OUTPUT, "Fourier_timeseries_H.pt"), weights_only=False).T
    min_len = min(C.shape[0], A.shape[0], H.shape[0])
    data = torch.stack([C[:min_len], A[:min_len], H[:min_len]], dim=0)
    times = torch.arange(min_len, dtype=torch.float32)
    return data, times


def load_protein_names():
    tsv_path = os.path.join(SANDBOX, "Step02_Clustering", "output", "downsampled_median.tsv")
    df = pd.read_csv(tsv_path, sep='\t', index_col=0)
    return list(df.index)


def run_model_inference(model, data, times, device):
    condition_names = ["Control", "H2O2", "Acetic_Acid"]
    OBS_LEN = 12
    latent_trajectories = {}
    pred_future = {}
    true_future = {}

    # PCA: encode FULL trajectory (all 16 points) → true latent, matching intepret.py Q2
    all_latent = []
    for i, name in enumerate(condition_names):
        full_x = data[i].unsqueeze(0).to(device)  # [1, T, N] all time points
        cond_id = torch.tensor([i], device=device)
        with torch.no_grad():
            lat = model.get_latent(full_x, cond_id, times.to(device))
        lat_np = lat.squeeze(0).cpu().numpy()
        latent_trajectories[name] = lat_np
        all_latent.append(lat_np)

    all_latent = np.concatenate(all_latent, axis=0)  # (3*T, latent_dim)
    pca = PCA(n_components=3)
    pca_components = pca.fit_transform(all_latent)

    # PC1 correlation with protein abundance
    protein_names = load_protein_names()
    all_abundance = data[:, :, :].reshape(-1, data.shape[2]).numpy()

    pc1 = pca_components[:, 0]
    correlations = []
    for p_idx in range(data.shape[2]):
        corr = np.corrcoef(pc1, all_abundance[:, p_idx])[0, 1]
        correlations.append(corr)

    pc1_df = pd.DataFrame({"gene_name": protein_names, "correlation_with_PC1": correlations})
    pc1_df.to_csv(os.path.join(OUTPUT, "gsea_ranked_list_PC1.csv"), index=False)

    # Extrapolation R²
    with torch.no_grad():
        for i, name in enumerate(condition_names):
            x_obs = data[i, :OBS_LEN, :].unsqueeze(0).to(device)
            cond_id = torch.tensor([i], device=device)
            pred, _ = model(x_obs, cond_id, times[OBS_LEN:].to(device))
            pred_future[name] = pred.squeeze(0).cpu().numpy()
            true_future[name] = data[i, OBS_LEN:, :].numpy()


def compute_gradients(model, data, times, device):
    OBS_LEN = 12
    condition_names = ["Control", "Acetic", "H2O2"]
    gradients = {}

    for i, name in enumerate(condition_names):
        x_obs = data[i, :OBS_LEN, :].unsqueeze(0).to(device)
        cond_id = torch.tensor([i], device=device)

        # Full sequence input with gradients enabled
        full_input = data[i].unsqueeze(0).to(device).clone().detach().requires_grad_(True)

        model.train()  # enable gradient flow
        mu, logvar = model.encoder(full_input)
        z0 = mu  # deterministic for gradients
        ode_f = functools.partial(model.ode_func, cond_emb=model.dynamic_embedding(cond_id))
        z_t = odeint_adjoint(ode_f, z0, times.to(device),
                             method='dopri5', adjoint_params=tuple(model.ode_func.parameters()))
        y_pred = model.decoder(z_t.permute(1, 0, 2))
        y_pred = y_pred + model.shift_embedding(cond_id).unsqueeze(1)

        loss = torch.nn.functional.mse_loss(y_pred, full_input)
        loss.backward()
        grads = full_input.grad.squeeze(0).cpu().numpy()  # (T, Nprot)
        model.eval()

        df = pd.DataFrame(grads.T, columns=[f"T{t}" for t in range(grads.shape[0])])
        df.index = load_protein_names()
        gradients[name] = df

    # Differential gradients
    for pair, (a, b) in [("Acetic_vs_Control", ("Acetic", "Control")),
                           ("H2O2_vs_Control", ("H2O2", "Control"))]:
        common = sorted(set(gradients[a].index) & set(gradients[b].index))
        # Ensure same columns
        cols = [c for c in gradients[a].columns if c in gradients[b].columns]
        diff = gradients[a].loc[common, cols] - gradients[b].loc[common, cols]
        path = os.path.join(OUTPUT, f"Diff_Gradients_{pair}_temporal_gradients.csv")
        diff.to_csv(path)


def call_r_script():
    r_script = os.path.join(BASE, "Step05_GSEA.R")
    ret = subprocess.run(["Rscript", r_script], cwd=BASE, capture_output=True,
                         text=True, encoding='utf-8', errors='replace')
    print(ret.stdout)
    if ret.returncode != 0:
        print("R stderr:", ret.stderr)
        print("[WARN] R script had errors")


def main():
    print("Step05 — Interpretation...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    preprocess_data()
    data, times = load_tensors()

    model_path = os.path.join(SANDBOX, "Step04_Modeling", "model", "best_model.pth")
    ckpt = torch.load(model_path, map_location=device, weights_only=False)
    if 'model_params' in ckpt:
        params = {k: v for k, v in ckpt['model_params'].items()
                  if not k.startswith('data_')}
    else:
        params = {'input_dim': data.shape[2], 'latent_dim': 128, 'hidden_dim': 256,
                  'num_conditions': data.shape[0], 'condition_dim': 32}
    model = GAM_LatentODE(**params).to(device)
    model.load_state_dict(ckpt['model_state_dict'])
    model.reparameterize = lambda mu, logvar: mu
    model.eval()

    run_model_inference(model, data, times, device)
    compute_gradients(model, data, times, device)
    call_r_script()

    print("Step05 — Done.")


if __name__ == "__main__":
    main()
