# Step3_efig4.py — EFig.4d/e: GAM MaxDiff vs canonical ESR validation
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

FILE_GAM = "../output/GAM_maxdiff_ranking.csv"
FILE_ESR = "../../ESR_gene_parameters.xls"

os.makedirs("../output", exist_ok=True)
os.makedirs("../figure", exist_ok=True)

gam = pd.read_csv(FILE_GAM)
gam.columns = gam.columns.str.strip()
gam['ProteinID'] = gam['ProteinID'].astype(str).str.strip().str.upper()

esr = pd.read_excel(FILE_ESR)[['Name', 'ESR']].dropna(subset=['ESR'])
esr['Name'] = esr['Name'].astype(str).str.strip().str.upper()

merged = pd.merge(gam, esr, left_on='ProteinID', right_on='Name', how='left')
merged['ESR'] = merged['ESR'].fillna('Non-ESR')


def plot_esr_validation(df, metric_col, title_prefix, output_filename):
    ranked = df.sort_values(by=metric_col, ascending=False).reset_index(drop=True)
    ranked['Rank'] = ranked.index
    ranks_up   = ranked[ranked['ESR'] == 'up']['Rank'].values
    ranks_down = ranked[ranked['ESR'] == 'down']['Rank'].values

    fig, (ax_curve, ax_rug, ax_dens) = plt.subplots(
        3, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={'height_ratios': [2, 1, 1.5]})
    plt.subplots_adjust(hspace=0.05)

    # Panel A: MaxDiff curve
    x, y = ranked['Rank'], ranked[metric_col]
    ax_curve.plot(x, y, color='black', lw=1.5, zorder=10)
    ax_curve.fill_between(x, y, 0, where=(y > 0),  color='#D62728', alpha=0.3, linewidth=0)
    ax_curve.fill_between(x, y, 0, where=(y <= 0), color='#2CA02C', alpha=0.3, linewidth=0)
    ax_curve.set_ylabel("GAM MaxDiff\n(Treatment - Control)", fontsize=12, fontweight='bold')
    ax_curve.set_title(f"Validation: {title_prefix} MaxDiff vs. Canonical ESR",
                       fontsize=14, pad=15, fontweight='bold')
    ax_curve.axhline(0, color='black', linestyle='--', lw=0.8)
    ax_curve.text(0, y.max(), f' Stress Induced\n({title_prefix} > Control)',
                  ha='left', va='top', fontweight='bold', color='#8B0000', fontsize=11)
    ax_curve.text(len(x), y.min(), f'Growth Repressed\n({title_prefix} < Control)',
                  ha='right', va='bottom', fontweight='bold', color='#006400', fontsize=11)

    # Panel B: Rug
    ax_rug.vlines(ranks_up,   0.5, 1, color='#D62728', lw=0.6, alpha=0.8)
    ax_rug.vlines(ranks_down, 0, 0.5, color='#2CA02C', lw=0.6, alpha=0.8)
    ax_rug.set_ylim(0, 1); ax_rug.set_yticks([0.25, 0.75])
    ax_rug.set_yticklabels(['ESR\nRepressed', 'ESR\nInduced'], fontsize=10, fontweight='bold')
    for sp in ['bottom','top']: ax_rug.spines[sp].set_visible(False)

    # Panel C: Density
    sns.kdeplot(ranks_up,   ax=ax_dens, color='#D62728', fill=True, alpha=0.3, lw=2, label='Stress Genes')
    sns.kdeplot(ranks_down, ax=ax_dens, color='#2CA02C', fill=True, alpha=0.3, lw=2, label='Growth Genes')
    ax_dens.set_xlim(0, len(ranked))
    ax_dens.set_ylabel("Gene Density", fontsize=12, fontweight='bold')
    ax_dens.set_xlabel(f"Rank in {title_prefix} Response (Up → Down)", fontsize=12, fontweight='bold')
    ax_dens.legend(loc='upper center', frameon=False, fontsize=11)

    for ax in [ax_curve, ax_rug, ax_dens]:
        ax.grid(axis='x', linestyle=':', alpha=0.5, color='gray')
        for sp in ax.spines.values():
            sp.set_edgecolor('black'); sp.set_linewidth(1.0)

    plt.savefig(output_filename, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)

plot_esr_validation(merged, 'MaxDiff_A_vs_C', 'Acetic Acid', '../figure/ef4d.pdf')
plot_esr_validation(merged, 'MaxDiff_H_vs_C', 'H2O2', '../figure/ef4e.pdf')
