# Fig5b.py — PC1 ESR validation (curve + rug + density)

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os

BASE = os.path.dirname(os.path.abspath(__file__))

file_pc1 = os.path.join(BASE, "..", "output", "gsea_ranked_list_PC1.csv")
file_esr = os.path.join(BASE, "..", "..", "ESR_gene_parameters.xls")

my_data = pd.read_csv(file_pc1)
my_data.columns = my_data.columns.str.strip()
if 'gene_name' not in my_data.columns:
    my_data.columns = ['gene_name', 'correlation_with_PC1']

esr_data = pd.read_excel(file_esr)
esr_subset = esr_data[['Name', 'ESR']].dropna(subset=['ESR'])

my_data['gene_name'] = my_data['gene_name'].astype(str).str.strip().str.upper()
esr_subset['Name'] = esr_subset['Name'].astype(str).str.strip().str.upper()
merged_df = pd.merge(my_data, esr_subset, left_on='gene_name', right_on='Name', how='left')
merged_df['ESR'] = merged_df['ESR'].fillna('Non-ESR')

ranked_df = merged_df.sort_values(by='correlation_with_PC1', ascending=False).reset_index(drop=True)
ranked_df['Rank'] = ranked_df.index

ranks_down = ranked_df[ranked_df['ESR'] == 'down']['Rank'].values
ranks_up = ranked_df[ranked_df['ESR'] == 'up']['Rank'].values

plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42})

fig, (ax_curve, ax_rug, ax_dens) = plt.subplots(3, 1, figsize=(10, 8), sharex=True,
                                                gridspec_kw={'height_ratios': [2, 1, 1.5]})
plt.subplots_adjust(hspace=0.05)

x = ranked_df['Rank']
y = ranked_df['correlation_with_PC1']

ax_curve.plot(x, y, color='black', lw=1.5)
ax_curve.fill_between(x, y, 0, where=(y > 0), color='green', alpha=0.3, interpolate=True)
ax_curve.fill_between(x, y, 0, where=(y <= 0), color='red', alpha=0.3, interpolate=True)
ax_curve.set_ylabel("PC1 Correlation", fontsize=12, fontweight='bold')
ax_curve.set_title("Validation: PC1 Axis captures the ESR Growth-Stress Continuum",
                   fontsize=14, pad=15)
ax_curve.axhline(0, color='black', linestyle='--', lw=0.8)
ax_curve.text(0, y.max(), ' Growth End\n(PC1 High)', ha='left', va='top',
              fontweight='bold', color='green')
ax_curve.text(len(x), y.min(), 'Stress End \n(PC1 Low)', ha='right', va='bottom',
              fontweight='bold', color='red')

ax_rug.vlines(ranks_down, 0.5, 1, color='#2CA02C', lw=0.5, alpha=0.8,
              label='ESR Repressed (Growth)')
ax_rug.vlines(ranks_up, 0, 0.5, color='#D62728', lw=0.5, alpha=0.8,
              label='ESR Induced (Stress)')
ax_rug.set_ylim(0, 1)
ax_rug.set_yticks([0.25, 0.75])
ax_rug.set_yticklabels(['ESR\nInduced', 'ESR\nRepressed'], fontsize=10, fontweight='bold')
ax_rug.legend(loc='center right', frameon=True)
ax_rug.spines['bottom'].set_visible(False)
ax_rug.spines['top'].set_visible(False)
ax_rug.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

sns.kdeplot(ranks_down, ax=ax_dens, color='#2CA02C', fill=True, alpha=0.3,
            label='Growth Genes', lw=2)
sns.kdeplot(ranks_up, ax=ax_dens, color='#D62728', fill=True, alpha=0.3,
            label='Stress Genes', lw=2)
ax_dens.set_xlim(0, len(ranked_df))
ax_dens.set_ylabel("Gene Density", fontsize=12, fontweight='bold')
ax_dens.set_xlabel("Rank in PC1 (Growth → Death)", fontsize=12, fontweight='bold')
ax_dens.legend(loc='upper center', frameon=False)

for ax in [ax_curve, ax_rug, ax_dens]:
    ax.grid(axis='x', linestyle=':', alpha=0.5)

out = os.path.join(BASE, "..", "figure", "fig5b.pdf")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, bbox_inches='tight', dpi=300)
plt.close(fig)
print("  -> fig5b.pdf")
