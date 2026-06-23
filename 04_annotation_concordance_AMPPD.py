"""
04_annotation_concordance_AMPPD.py
Author: Dimitrios Kyriakis

Cross-dataset pseudobulk correlation: miBRAIN organoids vs AMP-PD
postmortem brain reference to validate cell type annotation.

Usage:
    python 04_annotation_concordance_AMPPD.py \
        --query-h5ad    /path/to/miBRAIN_annotated.h5ad \
        --ref-h5ad      /path/to/AMPPD_reference.h5ad \
        --outdir        /path/to/output \
        --figure-dir    /path/to/figures

Arguments:
    --query-h5ad    Annotated miBRAIN AnnData (output of step 03)
    --ref-h5ad      AMP-PD / postmortem brain reference AnnData
                    (T. Clarence et al., FreshMG freeze 3)
    --outdir        Directory for output CSVs
    --figure-dir    Directory for output figures (PDF + PNG)
"""

import argparse
import os
import sys
import random
import warnings
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils_correlation import *

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="miBRAIN vs AMP-PD pseudobulk annotation concordance")
parser.add_argument("--query-h5ad",  required=True,
                    help="Annotated miBRAIN AnnData (.h5ad) from step 03")
parser.add_argument("--ref-h5ad",    required=True,
                    help="AMP-PD postmortem brain reference AnnData (.h5ad)")
parser.add_argument("--outdir",      required=True,
                    help="Output directory for correlation CSVs")
parser.add_argument("--figure-dir",  default="./figures",
                    help="Output directory for figures [default: ./figures]")
args = parser.parse_args()

h5ad1      = args.query_h5ad
h5ad2      = args.ref_h5ad
outdir     = args.outdir.rstrip("/") + "/"
figure_dir = args.figure_dir.rstrip("/") + "/"

os.makedirs(outdir,     exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
# Configure logging to stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    stream=sys.stderr
)

seed = 20210224
random.seed(seed)
np.random.seed(seed)
sc.settings.seed = seed
sc.settings.n_jobs = 1

data1 = pg.read_input(h5ad1)
data2 = pg.read_input(h5ad2)

data2_sub = data2[:, ~data2.var['gene_name'].duplicated(keep='first')].copy()
# 1. Save the current index (Ensembl IDs) into a new column
data2_sub.var['ENSEMBL'] = data2_sub.var.index
# 2. Set 'gene_name' as the new index (drops it from columns by default)
data2_sub.var.set_index('gene_name', inplace=True)
# 3. (Optional) Rename the index back to 'featurekey' if required by downstream tools
data2_sub.var.index.name = 'featurekey'
data2_sub.var['gene_name'] = data2_sub.var.index

# 2. Identify the overlapping genes
common_genes = data1.var_names.intersection(data2_sub.var_names)
logging.info(f"Retaining {len(common_genes)} shared genes.")

# --- Cluster to cell type mapping ---
# --- Cluster to cell type mapping ---
cluster_map = {
    '0':  'Mural/stromal cells',
    '1':  'Mural/stromal cells',
    '2':  'Neurons',
    '3':  'Early neurons',
    '4':  'Early neurons',
    '5':  'Astrocytes',
    '6':  'Progenitors/OPC',
    '7':  'Endothelial cells',
    '8':  'Perivascular progenitors',
    '9':  'Astrocytes',
    '10': 'Microglia',
}

# --- Curated palette — perceptually distinct, colorblind-friendly ---
cell_type_palette = {
    'Mural/stromal cells':    '#E64B35',  # red
    'Neurons':                 '#4DBBD5',  # cyan
    'Early neurons':      '#00A087',  # teal
    # 'Early neurons':         '#3C5488',  # navy
    'Astrocytes':                         '#FF9900',  # amber
    'Progenitors/OPC':       '#8491B4',  # slate
    'Endothelial cells':                  '#91D1C2',  # mint
    'Perivascular progenitors':           '#B09C85',  # taupe
    'Immature astroglia':                 '#F7DC6F',  # yellow
    'Microglia':                          '#7E6148',  # brown
}


# 1. Map labels and fill missing values with the original column
data1.obs['cell_type'] = (
    data1.obs['leiden_res0_20']
    .astype(str)
    .map(cluster_map)
    .fillna(data1.obs['leiden_res0_20'].astype(str))
    .astype('category')
)

data1.obs['class'] = data1.obs['cell_type']

data2_sub = data2[:, ~data2.var['gene_name'].duplicated(keep='first')].copy()
# 1. Save the current index (Ensembl IDs) into a new column
data2_sub.var['ENSEMBL'] = data2_sub.var.index

# 2. Set 'gene_name' as the new index (drops it from columns by default)
data2_sub.var.set_index('gene_name', inplace=True)

# 3. (Optional) Rename the index back to 'featurekey' if required by downstream tools
data2_sub.var.index.name = 'featurekey'
data2_sub.var['gene_name'] = data2_sub.var.index

# 2. Identify the overlapping genes
common_genes = data1.var_names.intersection(data2_sub.var_names)
logging.info(f"Retaining {len(common_genes)} shared genes.")

# 3. Slice the objects securely
data1_common = data1[:, common_genes].copy()
data2_common = data2_sub[:, common_genes].copy()



# # Calculate aggregation

group_by1 = 'cell_type'
group_by2 = 'class'
mat_key1  = 'counts'
mat_key2  = 'X'
prefix    = 'AMPPD_class'
logging.info("Aggregating pseudobulk matrices...")
agr_mat_h5ad1 = pge.pb_agg_by_cluster_loaded(data1_common, group_by1, log1p=True, mat_key=mat_key1)
agr_mat_h5ad2 = pge.pb_agg_by_cluster_loaded(data2_common, group_by2, log1p=True, mat_key=mat_key2)

subset2   =  True
degs_calc = 'both' 
nominal   =  False
deg_intersection=True
min_variance_percentile=75
scale     = 'robust'
scale_before_spearman = False
p_thresh = 0.05

# # DEGS h5ad1

adata1 =  data1_common.to_anndata()
adata1.X = adata1.layers['counts'].copy()
sc.pp.normalize_total(adata1, target_sum=1e4)
sc.pp.log1p(adata1)

# 4. Compute DEGs (Wilcoxon)
sc.tl.rank_genes_groups(adata1, group_by1, method='wilcoxon',use_raw=False)

result = adata1.uns['rank_genes_groups']
groups = result['names'].dtype.names

all_sig_genes = set()
for group in groups:
    mask = result['pvals_adj'][group] < p_thresh  # FDR correction applied here
    all_sig_genes.update(result['names'][group][mask])
degs1 = list(all_sig_genes)

# # DEGS h5ad2

adata2 =  data2_common.to_anndata()

max_cells_per_group=3000
np.random.seed(seed)
counts = adata2.obs[group_by2].value_counts()
keep_clusters = counts[counts >= 2].index.tolist()
sampled_indices = []
for cluster in keep_clusters:
    cluster_indices = adata2.obs_names[adata2.obs[group_by2] == cluster]
    if len(cluster_indices) > max_cells_per_group:
        choice = np.random.choice(cluster_indices, max_cells_per_group, replace=False)
        sampled_indices.extend(choice)
    else:
        sampled_indices.extend(cluster_indices)
adata2 = adata2[sampled_indices].copy()

sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)

# 4. Compute DEGs (Wilcoxon)
sc.tl.rank_genes_groups(adata2, group_by2, method='wilcoxon',use_raw=False)

result = adata2.uns['rank_genes_groups']
groups = result['names'].dtype.names

all_sig_genes = set()
for group in groups:
    mask = result['pvals_adj'][group] < p_thresh  # FDR correction applied here
    all_sig_genes.update(result['names'][group][mask])
degs2 = list(all_sig_genes)

# 1. Create the master list (union) of all significant genes
union_degs = set(degs1).union(degs2)

# 2. Extract the available genes (columns) from both aggregated matrices
available_genes_1 = set(agr_mat_h5ad1.columns)
available_genes_2 = set(agr_mat_h5ad2.columns)

# 3. Intersect the DEG union with the genes physically present in BOTH matrices
valid_common_genes = sorted(
    union_degs.intersection(available_genes_1).intersection(available_genes_2)
)

mat1_b_scale_agr = agr_mat_h5ad1.loc[:, valid_common_genes]
mat2_b_scale_agr = agr_mat_h5ad2.loc[:, valid_common_genes]

mat1_a_scale_agr = (mat1_b_scale_agr - mat1_b_scale_agr.mean(axis=0)) / mat1_b_scale_agr.std(axis=0)
mat1_a_scale_agr.dropna(axis=1, how="any", inplace=True)

mat2_a_scale_agr = (mat2_b_scale_agr - mat2_b_scale_agr.mean(axis=0)) / mat2_b_scale_agr.std(axis=0)
mat2_a_scale_agr.dropna(axis=1, how="any", inplace=True)

corr_pearson = pge.pearson_corr(mat1_a_scale_agr, mat2_a_scale_agr)
corr_spearman = pge.spearman_corr(mat1_a_scale_agr, mat2_a_scale_agr)

plot_correlation(corr_spearman.T,
    prefix='AMPPD',
    method='spearman',
    fig_width=7,
    fig_height=5,
    save_plot=True,
    figure_path=figure_dir,
    annot=False)

plot_correlation(corr_pearson.T,
    prefix='AMPPD',
    method='pearson',
    fig_width=7,
    fig_height=5,
    save_plot=True,
    figure_path=figure_dir,
    annot=False)

corr_pearson.to_csv(outdir + 'AMPPD_Pearson_Correlation.csv')
corr_spearman.to_csv(outdir + 'AMPPD_Spearman_Correlation.csv')

cor_df_ordered = corr_spearman.T

cor_df_ordered


corr_pearson =  pd.read_csv('AMPPD_Pearson_Correlation.csv')
corr_spearman = pd.read_csv('AMPPD_Spearman_Correlation.csv',index_col='barcodekey')


output_pdf = './figures/AMPPD_Spearman_Correlation.pdf'
# Reorder rows so matching clusters fall on the diagonal
# Columns: Astro, EN, Endo, Mural, Myeloid, OPC
row_order = ['Astro', 'OPC', 'EN', 'Endo', 'Myeloid',  'Mural']
col_order  = ['Astrocytes', 'Progenitors/OPC','Early neurons','Neurons', 'Endothelial cells', 
               'Microglia','Perivascular progenitors','Mural/stromal cells']
df = cor_df_ordered.T
df_ordered = df.loc[row_order, col_order]


fig, ax = plt.subplots(figsize=(5, 5))

sns.heatmap(
    df_ordered, 
    annot=False, 
    cmap='RdBu_r',   # Red/Blue diverging colormap highlights positive/negative separation
    center=0,        # Centered at 0 for correlation-like data
    fmt=".2f",
    linewidths=0.5,
    linecolor  = "black",
    cbar_kws={'label': 'spearman r'},
    ax=ax
)

ax.set_title('Cell Type Annotation Concordance')
ax.set_ylabel('Human Postmortem Brain')#T.Clarence et.al.')
ax.set_xlabel('miBRAIN')
# Format tick orientation cleanly
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
plt.tight_layout()
plt.show()
fig.savefig(output_pdf, dpi=300)
plt.close(fig)
# logging.info("Pipeline step complete.")
df_ordered.T.to_csv('2026_06_AMPPD_Spearman_Correlation_sub.csv')
