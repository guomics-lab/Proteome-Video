# Step2_efig3.py — EFig.3c: ESR enrichment Fisher test and bar chart
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os, sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

FILE_MASIGPRO  = os.path.join(script_dir, "..", "output", "maSigPro_Final_Results", "maSigPro_Analysis_Results.xlsx")
FILE_WGCNA     = os.path.join(script_dir, "..", "output", "WGCNA_module_assignments.csv")
FILE_BACKGROUND = os.path.join(script_dir, "..", "..", "Step00_Preprocess", "output", "imputed_temporal_matrix.tsv")
FILE_ESR        = os.path.join(script_dir, "..", "..", "ESR_gene_parameters.xls")

# Ensure output and figure directories exist
os.makedirs(os.path.join(script_dir, "..", "output"), exist_ok=True)
os.makedirs(os.path.join(script_dir, "..", "figure"), exist_ok=True)

# ---- Load background universe ----
df_bg = pd.read_csv(FILE_BACKGROUND, sep='\t')
df_bg['Genes'] = df_bg['Genes'].astype(str).str.upper().str.strip()
universe_genes = set(df_bg['Genes'])
N = len(universe_genes)

# ---- Load ESR reference ----
df_esr = pd.read_excel(FILE_ESR)
df_esr = df_esr[['Name', 'ESR']].dropna(subset=['ESR'])
df_esr['Name'] = df_esr['Name'].astype(str).str.upper().str.strip()
raw_esr_up   = set(df_esr[df_esr['ESR'] == 'up']['Name'])
raw_esr_down = set(df_esr[df_esr['ESR'] == 'down']['Name'])
bg_esr_up   = universe_genes.intersection(raw_esr_up)
bg_esr_down = universe_genes.intersection(raw_esr_down)
K_up, K_down = len(bg_esr_up), len(bg_esr_down)


def run_fisher(gene_groups, name_col, group_col, output_csv):
    """Run Fisher exact test for each group against ESR up/down."""
    results = []
    for group_name, group_genes in gene_groups.items():
        n = len(group_genes)
        # vs ESR Up
        hit_up = len(group_genes.intersection(bg_esr_up))
        a_up, b_up = hit_up, n - hit_up
        c_up, d_up = K_up - hit_up, N - n - (K_up - hit_up)
        _, p_up = fisher_exact([[a_up, b_up], [c_up, d_up]], alternative='greater')
        # vs ESR Down
        hit_down = len(group_genes.intersection(bg_esr_down))
        a_down, b_down = hit_down, n - hit_down
        c_down, d_down = K_down - hit_down, N - n - (K_down - hit_down)
        _, p_down = fisher_exact([[a_down, b_down], [c_down, d_down]], alternative='greater')
        results.append({
            name_col: group_name, 'Size': n,
            'ESR_Up_Hits': hit_up, 'p_Up': p_up,
            'ESR_Down_Hits': hit_down, 'p_Down': p_down
        })
    df = pd.DataFrame(results)
    # FDR correction
    all_p = df['p_Up'].tolist() + df['p_Down'].tolist()
    valid = [p for p in all_p if pd.notna(p)]
    if valid:
        _, corrected, _, _ = multipletests(valid, method='fdr_bh')
        pmap = dict(zip(valid, corrected))
        df['adj_p_Up']   = df['p_Up'].map(pmap)
        df['adj_p_Down'] = df['p_Down'].map(pmap)
    df.to_csv(output_csv, index=False)
    return df


# ========== 1. maSigPro Fisher test ==========
df_mas = pd.read_excel(FILE_MASIGPRO)
df_mas.columns = df_mas.columns.str.strip()
df_mas['Gene'] = df_mas['Gene'].astype(str).str.upper().str.strip()
df_mas['Cluster'] = df_mas['Cluster'].astype(str).str.strip()

# Merge Cluster_2,4,6
merge_targets = ['Cluster_2', 'Cluster_4', 'Cluster_6']
merged_name = 'Cluster_2_4_6_Merged'

# Track raw count BEFORE merging (for bar chart n=591)
merged_raw_n = int(df_mas[df_mas['Cluster'].isin(merge_targets)].shape[0])

def assign_group(c):
    return merged_name if c in merge_targets else c
df_mas['Group'] = df_mas['Cluster'].apply(assign_group)

mas_groups = {g: set(df_mas[df_mas['Group'] == g]['Gene'])
              for g in sorted(df_mas['Group'].unique())}
mas_df = run_fisher(mas_groups, 'Group', 'Group',
                    os.path.join(script_dir, "..", "output", 'ef3c_maSigPro_SourceData.csv'))

# ========== 2. WGCNA Fisher test ==========
df_wgcna = pd.read_csv(FILE_WGCNA)
df_wgcna['Gene'] = df_wgcna['Gene'].astype(str).str.upper().str.strip()
df_wgcna = df_wgcna.assign(Gene=df_wgcna['Gene'].str.split(';')).explode('Gene')
df_wgcna = df_wgcna.drop_duplicates().reset_index(drop=True)
df_wgcna['ModuleColor'] = df_wgcna['ModuleColor'].astype(str).str.strip()

wgcna_groups = {m: set(df_wgcna[df_wgcna['ModuleColor'] == m]['Gene'])
                for m in sorted(df_wgcna['ModuleColor'].unique())}
wgcna_df = run_fisher(wgcna_groups, 'Module', 'ModuleColor',
                      os.path.join(script_dir, "..", "output", 'ef3c_WGCNA_SourceData.csv'))

# ========== 3. Bar chart (ef3c.pdf) ==========
background_ratio = (K_down / N) * 100

# Read computed values dynamically
data_rows = []
# maSigPro merged cluster (n from raw pre-merge count)
row_mas = mas_df[mas_df['Group'] == merged_name].iloc[0]
data_rows.append({
    "name": "Cluster 2,4,6\n(Combined)",
    "n": merged_raw_n, "k": int(row_mas['ESR_Down_Hits']),
    "p": row_mas['p_Down']
})
# Red module
row_red = wgcna_df[wgcna_df['Module'] == 'red'].iloc[0]
data_rows.append({
    "name": "Red Module",
    "n": int(row_red['Size']), "k": int(row_red['ESR_Down_Hits']),
    "p": row_red['p_Down']
})
# Turquoise module
row_turq = wgcna_df[wgcna_df['Module'] == 'turquoise'].iloc[0]
data_rows.append({
    "name": "Turquoise Module",
    "n": int(row_turq['Size']), "k": int(row_turq['ESR_Down_Hits']),
    "p": row_turq['p_Down']
})

df = pd.DataFrame(data_rows)
df['Observed_Ratio'] = (df['k'] / df['n']) * 100
df['LogP'] = -np.log10(df['p'].replace(0, np.nan))

df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['Observed_Ratio', 'LogP'])
if len(df) == 0:
    print("ERROR: No valid data rows for bar chart"); sys.exit(1)
df = df.reset_index(drop=True)

fig, ax = plt.subplots(figsize=(6, 3.5))
norm = mcolors.Normalize(vmin=1.3, vmax=8)
cmap = plt.cm.Reds

bars = ax.barh(df.index, df['Observed_Ratio'], color=cmap(norm(df['LogP'])),
               edgecolor='black', height=0.6, zorder=3)
ax.axvline(x=background_ratio, color='gray', linestyle='--', linewidth=2, zorder=2)
ax.text(background_ratio, -0.6, f'Random Expectation\n({background_ratio:.1f}%)',
        color='gray', fontsize=8, ha='center', va='top')

max_width = df['Observed_Ratio'].max()
for i, bar in enumerate(bars):
    width = bar.get_width()
    row = df.iloc[i]
    count_text = f"{int(row['k'])}/{int(row['n'])}"
    text_color = 'white' if row['LogP'] > 2 else 'black'
    if width > 5:
        ax.text(width / 2, bar.get_y() + bar.get_height() / 2, count_text,
                ha='center', va='center', color=text_color, fontweight='bold', fontsize=10)
    else:
        ax.text(width + 1, bar.get_y() + bar.get_height() / 2, count_text,
                ha='left', va='center', color='black', fontsize=10)
    p_str = f"$P$={row['p']:.1e}" if row['p'] < 0.001 else f"$P$={row['p']:.3f}"
    ax.text(width + 2, bar.get_y() + bar.get_height() / 2, p_str,
            ha='left', va='center', color='black', fontsize=10)

ax.set_yticks(df.index)
ax.set_yticklabels(df['name'], fontsize=11, fontweight='bold', color='black')
ax.invert_yaxis()
ax.set_xlabel("Percentage of ESR Down Genes (%)", fontsize=10)
ax.set_xlim(0, max_width * 1.45)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.18, aspect=40, shrink=0.6)
cbar.set_label("-log10(P-value) Significance", fontsize=9)
cbar.ax.tick_params(labelsize=8)
plt.title("Fisher's Exact Test Results\n(Observed vs. Expected)", fontweight='bold', pad=10)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "..", "figure", "ef3c.pdf"), bbox_inches='tight')
print("  -> ef3c.pdf")
