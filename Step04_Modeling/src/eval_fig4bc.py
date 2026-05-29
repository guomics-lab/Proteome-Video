# eval_fig4bc.py — Fig.4b/c: Latent ODE inference and visualization

import torch, torch.nn as nn, numpy as np, os, sys, json, functools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
from torchdiffeq import odeint_adjoint
import seaborn as sns

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
        ce = cond_emb.expand(z.shape[0], -1)
        return self.baseline_net(z) + self.diff_net(torch.cat([z, ce], dim=1))

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
        return mu + 0.5 * logvar.exp().sqrt() * torch.randn_like(mu)
    def forward(self, x_obs, cond_id, time_to_predict):
        mu, logvar = self.encoder(x_obs)
        z0 = self.reparameterize(mu, logvar)
        kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        dyn = self.dynamic_embedding(cond_id)
        ode_f = functools.partial(self.ode_func, cond_emb=dyn)
        z_t = odeint_adjoint(ode_f, z0, time_to_predict.to(x_obs.device),
                             method='dopri5', adjoint_params=tuple(self.ode_func.parameters()))
        z_traj = z_t.permute(1, 0, 2)
        y = self.decoder(z_traj) + self.shift_embedding(cond_id).unsqueeze(1)
        return y, kl
    def get_latent_trajectory(self, x_obs, cond_id, time_points):
        mu, _ = self.encoder(x_obs); z0 = mu
        dyn = self.dynamic_embedding(cond_id)
        ode_f = functools.partial(self.ode_func, cond_emb=dyn)
        z_t = odeint_adjoint(ode_f, z0, time_points, method='dopri5',
                             adjoint_params=tuple(self.ode_func.parameters()))
        return z_t.permute(1, 0, 2)

# ── Data ────────────────────────────────────────────────────────────────────
def load_raw_data(data_dir="data"):
    C = torch.load(os.path.join(data_dir, "Fourier_timeseries_C.pt"), weights_only=False).T
    A = torch.load(os.path.join(data_dir, "Fourier_timeseries_A.pt"), weights_only=False).T
    H = torch.load(os.path.join(data_dir, "Fourier_timeseries_H.pt"), weights_only=False).T
    s = [C, A, H]; m = min(x.shape[0] for x in s)
    full = torch.stack([x[:m, :] for x in s], dim=0)
    return full, torch.arange(m, dtype=torch.float32)

# ── Fig4b: per-condition scatter ────────────────────────────────────────────
def plot_scatter(y_true, y_pred, cond_name, output_path, color):
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42})
    yt, yp = y_true.flatten(), y_pred.flatten()
    r2, mse = r2_score(yt, yp), float(np.mean((yt - yp) ** 2))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(yt, yp, c=color, alpha=0.3, s=15, edgecolors='none', rasterized=True)
    lo, hi = min(yt.min(), yp.min()), max(yt.max(), yp.max())
    ax.plot([lo, hi], [lo, hi], 'k--', lw=1.5, alpha=0.8, label='Perfect Prediction')
    ax.set_xlabel('True Values', fontsize=16); ax.set_ylabel('Predicted Values', fontsize=16)
    ax.set_title(f'{cond_name}\n$R^2 = {r2:.4f}, MSE = {mse:.5f}$', fontsize=18)
    ax.legend(loc='upper left', frameon=True, fontsize=12)
    ax.grid(True, linestyle='-', alpha=0.3, color='lightgray')
    ax.set_aspect('equal', adjustable='box')
    fig.savefig(output_path, bbox_inches='tight', dpi=300); plt.close()
    return r2, mse

# ── Fig4c: PCA space visualisation ──────────────────────────────────────────
def plot_space_pca(data_per_condition, space_name, output_dir, colors_dict):
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 12,
                         'font.family': 'sans-serif', 'axes.linewidth': 1.2})
    names = list(data_per_condition.keys())
    all_data = np.concatenate([data_per_condition[n] for n in names], axis=0)
    emb = PCA(n_components=2, random_state=42).fit_transform(all_data)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.grid(False)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    si = 0
    for nm in names:
        n_pts = len(data_per_condition[nm]); pts = emb[si:si + n_pts]; si += n_pts
        if 'Control' in nm:  c, lbl = colors_dict['C'], 'Control'
        elif 'H2O2' in nm:   c, lbl = colors_dict['H'], 'H2O2'
        else:                c, lbl = colors_dict['A'], 'Acetic Acid'
        ax.plot(pts[:, 0], pts[:, 1], color=c, alpha=0.4, lw=1.5, zorder=1)
        ax.scatter(pts[:, 0], pts[:, 1], s=40, color=c, alpha=1.0, edgecolors='none', zorder=2, label=lbl)
        for t in [0, 3, 6, 9, 12, 15]:
            if t < n_pts:
                ax.annotate(str(t), (pts[t, 0], pts[t, 1]), textcoords='offset points',
                            xytext=(3, 3), color='black', fontsize=8, zorder=3)
    ax.set_xlabel('PCA Dim 1', fontweight='bold'); ax.set_ylabel('PCA Dim 2', fontweight='bold')
    h, l = ax.get_legend_handles_labels(); ax.legend(dict(zip(l, h)).values(), dict(zip(l, h)).keys(), frameon=False, fontsize=10)
    plt.tight_layout()
    fp = os.path.join(output_dir, f'{space_name}_PCA.pdf')
    fig.savefig(fp, bbox_inches='tight', dpi=300); plt.close()

# ── Main ────────────────────────────────────────────────────────────────────
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__)); os.chdir(script_dir)
    MODEL_PATH = os.path.join(script_dir, "..", "model", "best_model.pth")
    DATA_DIR, OBS_LEN = "../output", 12
    EVAL_DIR = os.path.join(script_dir, "..", "figure")
    os.makedirs(EVAL_DIR, exist_ok=True)

    if not os.path.exists(MODEL_PATH): print(f"ERROR: {MODEL_PATH}"); sys.exit(1)
    if not os.path.exists(DATA_DIR): print(f"ERROR: {DATA_DIR}"); sys.exit(1)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    colors = {"C": "#2E8B57", "H": "#1E90FF", "A": "#FF6347"}
    cond_names = ["Control", "H2O2", "Acetic_Acid"]

    full_data, full_time = load_raw_data(DATA_DIR)

    ckpt = torch.load(MODEL_PATH, map_location=device)
    if 'model_params' in ckpt:
        params = {k: v for k, v in ckpt['model_params'].items()
                  if k in ('input_dim','latent_dim','hidden_dim','num_conditions','condition_dim')}
    else:
        params = {'input_dim': int(full_data.shape[2]), 'latent_dim': 128,
                  'hidden_dim': 256, 'num_conditions': 3, 'condition_dim': 32}

    model = GAM_LatentODE(**params).to(device)
    model.load_state_dict(ckpt['model_state_dict'])
    model.reparameterize = lambda mu, logvar: mu  # deterministic
    model.eval()

    observed = full_data[:, :OBS_LEN, :]

    time_future = full_time[OBS_LEN:]

    # ── Fig4b: per-condition scatter ──
    with torch.no_grad():
        for i, nm in enumerate(cond_names):
            x_obs = observed[i].unsqueeze(0).to(device)
            cid = torch.tensor([i], device=device)
            y_pred_future, _ = model(x_obs, cid, time_future.to(device))
            y_pred_future = y_pred_future.squeeze(0).cpu().numpy()
            y_true_future = full_data[i, OBS_LEN:, :].numpy()
            plot_scatter(y_true_future, y_pred_future, nm,
                         os.path.join(EVAL_DIR, f'fig4b_{nm}_simple_scatter.pdf'), list(colors.values())[i])

    # ── Fig4c: Protein_Space_PCA + Latent_Space_PCA ──
    protein_data, latent_data = {}, {}
    with torch.no_grad():
        for i, nm in enumerate(cond_names):
            x_obs = observed[i].unsqueeze(0).to(device)
            cid = torch.tensor([i], device=device)
            # Protein space: full true trajectory (matching original eval_fige5.py)
            protein_data[nm] = full_data[i, :, :].numpy()
            # Latent trajectory through full time
            lat = model.get_latent_trajectory(x_obs, cid, full_time.to(device))
            latent_data[nm] = lat.squeeze(0).cpu().numpy()

    plot_space_pca(protein_data, "fig4c_Protein_Space", EVAL_DIR, colors)
    plot_space_pca(latent_data, "fig4c_Latent_Space", EVAL_DIR, colors)
    print("\nDone.")

if __name__ == "__main__":
    main()
