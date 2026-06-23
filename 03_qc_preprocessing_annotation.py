"""
03_qc_preprocessing_annotation.py
Author: Dimitrios Kyriakis

QC, normalisation, dimensionality reduction, and cell type annotation
of iPSC-derived miniature brain organoid snRNA-seq data (miBRAIN).

Usage:
    python 03_qc_preprocessing_annotation.py \
        --cellranger-dir  /path/to/01_cellranger_counts \
        --outdir          /path/to/output_qc \
        --ensembl-csv     /path/to/human_ensembl_metadata.csv \
        --mitocarta-csv   /path/to/Human.MitoCarta3.0.csv \
        --ref-agg-dir     /path/to/reference_aggregations \
        --annot-h5ad      /path/to/final_annotated.h5ad

Arguments:
    --cellranger-dir   Root directory of Cell Ranger count outputs
                       (one subdirectory per sample containing outs/)
    --outdir           Directory for QC and intermediate h5ad files
    --ensembl-csv      Human Ensembl metadata CSV (BioMart export):
                         required columns: Gene name, Gene stable ID,
                         Chromosome/scaffold name, Gene type
                         Download: https://www.ensembl.org/biomart
    --mitocarta-csv    Human MitoCarta 3.0 CSV
                         Download: https://www.broadinstitute.org/mitocarta
    --ref-agg-dir      Directory containing reference pseudobulk aggregation
                         CSVs used for Leiden cluster annotation
    --annot-h5ad       Output path for final annotated AnnData (.h5ad)
"""

import argparse
import os
import sys

from utils import *
import session_info

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="miBRAIN snRNA-seq QC and annotation")
parser.add_argument("--cellranger-dir", required=True,
                    help="Root directory of Cell Ranger count outputs")
parser.add_argument("--outdir",         required=True,
                    help="Output directory for h5ad files")
parser.add_argument("--ensembl-csv",    required=True,
                    help="Path to human Ensembl metadata CSV (BioMart)")
parser.add_argument("--mitocarta-csv",  required=True,
                    help="Path to Human.MitoCarta3.0.csv")
parser.add_argument("--ref-agg-dir",    required=True,
                    help="Directory of reference pseudobulk aggregation CSVs")
parser.add_argument("--annot-h5ad",     required=True,
                    help="Output path for final annotated h5ad file")
args = parser.parse_args()

inpath      = args.cellranger_dir.rstrip("/") + "/"
outpath     = args.outdir.rstrip("/") + "/"
ensembl_csv = args.ensembl_csv
mito_csv    = args.mitocarta_csv
pge_agr_path = args.ref_agg_dir.rstrip("/") + "/"
annot_h5ad  = args.annot_h5ad

os.makedirs(outpath, exist_ok=True)

# ---------------------------------------------------------------------------
# Session info
# ---------------------------------------------------------------------------
import session_info
session_info.show()

dirs = [d for d in os.listdir(inpath) if os.path.isdir(os.path.join(inpath, d))]

dirs

for sample in dirs:
    read_10x_pegasus(inpath,outpath,sample)

import os
h5_files = [
    os.path.abspath(os.path.join(outpath, f))
    for f in os.listdir(outpath)
    if f.endswith(("-a.h5ad", "-b.h5ad")) and os.path.isfile(os.path.join(outpath, f))
]
print(h5_files)

# # Concatenate Batches

h5ad_output = outpath+'2026_New.h5ad'

adatas = [sc.read_h5ad(fp) for fp in h5_files]
import anndata
batch_names = [path.split("/")[-1].replace(".h5ad", "") for path in h5_files]
batch_names
adata = anndata.concat(adatas, join='outer', label='Sample', keys=batch_names, index_unique='-')
adata.write_h5ad(h5ad_output)

# # Add Metadata

h5ad_input = outpath+'2026_New.h5ad'
h5ad_output=h5ad_input.replace('.h5ad', '.meatadata.h5ad')
prefix=h5ad_input.replace('.h5ad', '')
add_metadata(input=h5ad_input,
             output=h5ad_output,
             prefix=prefix,
             ensembl_metadata_file=ensembl_csv,
             mitocarta_file=mito_csv)

# # QC

h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.qc.h5ad')
prefix = h5ad_input.replace('.h5ad', '')
data = pg.read_input(h5ad_input, file_type='h5ad', genome='GRCh38', modality='rna')
qc(input=h5ad_input, output=h5ad_output, prefix=prefix, mitocarta_file=mito_csv)

# # log1p_norm

h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.qc_norm.h5ad')
prefix = h5ad_input.replace('.h5ad', '')
log1p_norm(input=h5ad_input, output=h5ad_output)

# # sig_score

h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.qc_norm_sign.h5ad')
sig_score(input=h5ad_input, output=h5ad_output)

# # HVG

h5ad_input = h5ad_output
hvg = pge.scanpy_hvf_h5ad(h5ad_file=h5ad_input, flavor='cell_ranger', batch_key='Sample', 
                          n_top_genes=5000, min_mean=0.0125, max_mean=3, min_disp=0.5, robust_protein_coding=True)

# # First Pass

h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.qc_norm_sign_harmonized.h5ad')

res=1.2 
class_label='leiden_labels_res'+str(int(res*10))
n_pcs= 50 
batch='Sample' 
n_neighbors= 15 
tsne= False


data = pg.read_input(h5ad_input)
### set scanpy hvg as hvf
data.var.highly_variable_features = False
data.var.loc[data.var.index.isin(hvg),'highly_variable_features'] = True
### final value counts
print(data.var.highly_variable_features.value_counts())


### pca/harmony/umap/tsne
pg.pca(data, n_components=n_pcs)
pg.elbowplot(data)
npc = min(data.uns["pca_ncomps"], n_pcs)
print('Using %i components for PCA' % npc)
pg.regress_out(data, attrs=['n_counts','percent_mito','cycle_diff'])
pg.run_harmony(data, batch=batch, rep='pca_regressed', max_iter_harmony=20, n_comps=npc)
pg.neighbors(data, rep='pca_regressed_harmony', use_cache=False, K=100, n_comps=npc)
pg.leiden(data, rep='pca_regressed_harmony', resolution=res, class_label=class_label)
pg.umap(data, rep='pca_regressed_harmony', n_neighbors=n_neighbors, rep_ncomps=npc)
# figure
sc.pl.umap(data.to_anndata(), color=[class_label], legend_loc='on data', frameon=False, legend_fontsize=5, legend_fontoutline=1, title=[class_label], size=1, wspace=0, ncols=1,)
pg.write_output(data, h5ad_output, file_type='h5ad')

# # Find doublets

h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.qc_norm_harm_DBL_filt.h5ad')

data=pg.read_input(h5ad_input)

### find doublets
class_label='leiden_labels_res'+str(int(res*10))
pg.infer_doublets(data, channel_attr = 'Channel', clust_attr = class_label, plot_hist=None)
pg.mark_doublets(data)
pg.scatter(data, attrs='demux_type', basis='umap', dpi=150, return_fig=True).savefig('figures/doublets.png')
print(data.uns['pred_dbl_cluster'])
### doublet counts
dc = data.obs['demux_type'].value_counts().reset_index()
print(dc)
pct_dbl = dc.loc[dc['demux_type']=='doublet','count'] / np.sum(dc.loc[:,'count']) * 100
print('Doublets: %.2f%%' % pct_dbl)
### filter doublets
pg.qc_metrics(data, select_singlets=True)
pg.filter_data(data)

pg.write_output(data, h5ad_output, file_type='h5ad')

# # run_global_clustering_umap

# Path constructed automatically from outpath
h5ad_input = outpath + '2026_New.meatadata.qc.qc_norm.qc_norm_sign.qc_norm_sign_harmonized.qc_norm_harm_DBL_filt.h5ad'

# h5ad_input = h5ad_output
h5ad_output = h5ad_input.replace('.h5ad', '.run_global_clustering_umap.h5ad')
batch_key ='Sample'
run_global_clustering_umap(
        input_h5ad=h5ad_input,
        output_h5ad=h5ad_output,
        batch_key=batch_key,
        n_top_genes=6000,min_mean=0.0125,max_mean= 3.0,
    min_disp= 0.5,
        n_pcs=30,
        k_graph=100,
        n_neighbors_umap=100
    )

# # Annotate

adata = sc.read_h5ad(h5ad_output)

sample_map = {
    'Blanchard-wt-snRNAseq-a':      'wt',
    'Blanchard-wt-snRNAseq-b':      'wt',
    'Blanchard-A53T-E4-snRNAseq-a': 'A53T-E4',
    'Blanchard-A53T-E4-snRNAseq-b': 'A53T-E4',
    'Blanchard-A53T-E3-snRNAseq-a': 'A53T-E3',
    'Blanchard-A53T-E3-snRNAseq-b': 'A53T-E3'
}

# 2. Create the new column 'Condition'
adata.obs['Condition'] = adata.obs['Sample'].map(sample_map)

querys = pge.pb_agg_by_cluster(h5ad_file=h5ad_output,
                           cluster_label='leiden_res0_20',
                           robust_var_label=False,
                           log1p=True,mat_key='counts')
Mydata = querys

cell_agr_path = pge_agr_path  # set via --ref-agg-dir argument
agr_to_compare ={
    '2025_HgOrgAtlas_pbagg':'2025_HgOrgAtlas_pbagg.csv',
    '2025_HgOrgAtlas_hvf':'2025_HgOrgAtlas_pbagg_hvf.csv',
    'snMulti_dh_Subclass':'snMultiome_data_of_the_developing_human_neocortex_Subclass_pagg.csv',
    'AMPPD_Class':'AMPPD_freeze2_7_cx_rareMut_2024_12_16nolayers_class.csv',
    'AMPPD_Subclass':'AMPPD_freeze2_7_cx_rareMut_2024_12_16nolayers_subclass.csv',
    'PsychAD_Subclass':'2025_PsychAD_pbagg_subclass.csv'
}

for agr in agr_to_compare.keys(): 
    compute_cor(Mydata=Mydata,agr_name=agr,cluster="leiden_res0_20",agr_file=agr_to_compare[agr],pge_agr_path=pge_agr_path)

import logging

# Define your mapping: { 'cluster_id': 'Cell Type' }
cluster_map = {
    '0': 'MC',
    '1': 'MC',
    '2': 'Neurons',
    '3': 'NPC',
    '4': 'IP',
    '5': 'Astrocytes',
    '6': 'Neuroepithelium',
    '7': 'Endothelial',
    '8': 'PC',
    '9': 'IPC-Glia',
    '10': 'Microglia'
}

# 1. Map labels and fill missing values with the original column
adata.obs['cell_type'] = (
    adata.obs['leiden_res0_20']
    .map(cluster_map)
    .fillna(adata.obs['leiden_res0_20'])
    .astype('category')
)

# 2. Verify and log
if adata.obs['cell_type'].isnull().any():
    missing = adata.obs.loc[adata.obs['cell_type'].isnull(), 'leiden_res0_20'].unique()
    logging.warning(f"Clusters remains NaN even after fallback: {missing}")

# palette = sns.color_palette("colorblind", n_colors=num_categories).as_hex()
sc.pl.umap(adata, color='cell_type', 
           legend_loc='on data', frameon=False, 
           legend_fontsize=5, legend_fontoutline=1, size=1, 
           wspace=0, ncols=1,palette='Set1')

h5ad_output

cols_to_remove = [
        # '2025_HgOrgAtlas_pbagg_corr_annot_pnas', '2025_HgOrgAtlas_pbagg_corr_score_pnas',
        '2025_HgOrgAtlas_hvf_corr_annot_pnas', '2025_HgOrgAtlas_hvf_corr_score_pnas',
        'snMulti_dh_Subclass_corr_annot_pnas', 'snMulti_dh_Subclass_corr_score_pnas',
        'AMPPD_Class_corr_annot_pnas', 'AMPPD_Class_corr_score_pnas',
        # 'AMPPD_Subclass_corr_annot_pnas', 'AMPPD_Subclass_corr_score_pnas',
        'PsychAD_Subclass_corr_annot_pnas', 'PsychAD_Subclass_corr_score_pnas'
    ]
    
# Φιλτράρισμα μόνο των στηλών που υπάρχουν ήδη για αποφυγή σφαλμάτων
existing_cols = [c for c in cols_to_remove if c in adata.obs.columns]

adata.obs.drop(columns=existing_cols, inplace=True)

h5ad_output = annot_h5ad  # set via --annot-h5ad argument

adata.write(h5ad_output)

querys = pge.pb_agg_by_cluster(h5ad_file=h5ad_output,
                           cluster_label='cell_type',
                           robust_var_label=False,
                           log1p=True,mat_key='counts')
Mydata = querys

for agr in agr_to_compare.keys(): 
    compute_cor(Mydata=Mydata,agr_name=agr,cluster="cell_type",agr_file=agr_to_compare[agr],pge_agr_path=pge_agr_path)

def compute_cor(Mydata,agr_name,agr_file,pge_agr_path,cluster="manual_ct", figsize=(12, 5),legend_on=False):
    comp_agr_df = pd.read_csv(pge_agr_path+agr_file, index_col="barcodekey")

    # Remove zero-variance columns first
    Mydata = Mydata.loc[:, Mydata.std(axis=0) > 0]
    comp_agr_df = comp_agr_df.loc[:, comp_agr_df.std(axis=0) > 0]
    
    # Scale safely
    Mydata_scaled = (Mydata - Mydata.mean(axis=0)) / Mydata.std(axis=0)
    comp_agr_df_scaled = (comp_agr_df - np.mean(comp_agr_df, axis=0)) / np.std(comp_agr_df, axis=0)
    
    
    
    # Clean NaN/Inf after scaling
    Mydata_scaled = Mydata_scaled.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    # Do the same for comp_agr_df_scaled if needed
    comp_agr_df_scaled = comp_agr_df_scaled.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    
    # Sort by index if needed
    Mydata_scaled = Mydata_scaled.sort_index(axis=0)
    comp_agr_df_scaled = comp_agr_df_scaled.sort_index(axis=0)
    
    # Now run correlation
    corr_mat_pge = pge.pearson_corr(Mydata_scaled, comp_agr_df_scaled)
    # If corr_mat_pge has clusters as columns and ref annotations as rows:
    best_hits = corr_mat_pge.idxmax(axis=0)  # gives the ref annotation with highest correlation for each cluster
    best_scores = corr_mat_pge.max(axis=0)   # highest correlation score
    
    # Combine into a dataframe
    best_hits_df = pd.DataFrame({
        "best_annotation": best_hits,
        "best_score": best_scores
    })
    print(best_hits_df)

    # Make sure best_hits_df index is string to match adata.obs["leiden"]
    best_hits_df.index = best_hits_df.index.astype(str)
    
    # Now map annotations and scores
    adata.obs[agr_name+"_corr_annot_pnas"] = adata.obs[cluster].astype(str).map(best_hits_df["best_annotation"])
    adata.obs[agr_name+"_corr_score_pnas"] = adata.obs[cluster].astype(str).map(best_hits_df["best_score"])


    # Calculate dynamic height based on number of rows
    # min_height ensures small matrices don't look squashed
    num_rows = corr_mat_pge.shape[0]
    dynamic_height = max(5, num_rows * 0.5) 
    
    # Calculate dynamic width based on clusters (columns)
    num_cols = corr_mat_pge.shape[1]
    dynamic_width = max(8, num_cols * 0.8)

    # Now plot with adjusted figsize
    g = sns.clustermap(
        corr_mat_pge, 
        metric="correlation", 
        method="average",
        z_score=None, 
        cmap="vlag", 
        center=0, 
        figsize=(dynamic_width, dynamic_height)
    )
    g.fig.suptitle(agr_name+'_corr',fontsize=16, y=1.02)
    # g.ax_heatmap.set_title(agr_name+'_corr')


    # Optional: Rotate row labels if they are long
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
    plt.show()

cell_agr_path = pge_agr_path  # set via --ref-agg-dir argument
agr_to_compare ={
    '2025_HgOrgAtlas_hvf':'2025_HgOrgAtlas_pbagg_hvf.csv',
    'snMulti_dh_Subclass':'snMultiome_data_of_the_developing_human_neocortex_Subclass_pagg.csv',
    'AMPPD_Class':'AMPPD_freeze2_7_cx_rareMut_2024_12_16nolayers_class.csv',
    'PsychAD_Class': '2025_PsychAD_pbagg_hvf.csv'
}

for agr in agr_to_compare.keys(): 
    compute_cor(Mydata=Mydata,agr_name=agr,cluster="cell_type",agr_file=agr_to_compare[agr],pge_agr_path=pge_agr_path)

# # Percentages

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sample_map = {
    'Blanchard-wt-snRNAseq-a':      'Control',
    'Blanchard-wt-snRNAseq-b':      'Control',
    'Blanchard-A53T-E4-snRNAseq-a': 'A53T-E4',
    'Blanchard-A53T-E4-snRNAseq-b': 'A53T-E4',
    'Blanchard-A53T-E3-snRNAseq-a': 'A53T-E3',
    'Blanchard-A53T-E3-snRNAseq-b': 'A53T-E3'
}

# 2. Create the new column 'Condition'
adata.obs['Condition'] = adata.obs['Sample'].map(sample_map)

# 1. Set the categorical order: Control first, then mutants
# This ensures they appear in this specific order in the legend and on the bar clusters
categories = ['Control', 'A53T-E3', 'A53T-E4']
adata.obs['Condition'] = pd.Categorical(adata.obs['Condition'], categories=categories, ordered=True)

# 2. Calculate counts and percentages
# Normalized by the total size of each Condition (Columns sum to 100%)
counts = adata.obs.groupby(['cell_type', 'Condition'], observed=True).size().unstack(fill_value=0)
pcts = (counts / counts.sum() * 100).melt(ignore_index=False).reset_index()

# 3. Determine the sort order
# Sort by the abundance in the 'wt' condition (or total sum if preferred)
# Here I sort by total abundance across all samples so big clusters are first
cell_type_order = pcts.groupby('cell_type')['value'].sum().sort_values(ascending=False).index

# 4. Plot
plt.figure(figsize=(10, 4))
sns.barplot(
    data=pcts, 
    x='cell_type', 
    y='value', 
    hue='Condition', 
    order=cell_type_order, 
    palette={
        'Control': 'black',       # Control
        'A53T-E3': 'orange', # Mutant 1
        'A53T-E4': 'firebrick' # Mutant 2
    },
    edgecolor='black',       # Adds a nice border to bars
    linewidth=0.5
)

plt.ylabel('Percentage of Condition (%)')
plt.xlabel('cell_type')
plt.title('Cell Type Composition per Genotype')
plt.xticks(rotation=90)
plt.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left') # Move legend out
plt.tight_layout()
plt.show()

adata.obs['Condition'] = adata.obs['Sample'].map(sample_map)

# Ορισμός κατηγοριών εξαιρώντας το E4
target_categories = ['Control', 'A53T-E3']
logging.info(f"Filtering dataset to include: {target_categories}")

# Δημιουργία subset του obs για το γράφημα
plot_df = adata.obs[adata.obs['Condition'].isin(target_categories)].copy()
plot_df['Condition'] = pd.Categorical(plot_df['Condition'], categories=target_categories, ordered=True)

# 2. Calculate counts and percentages
# Χρησιμοποιούμε το φιλτραρισμένο dataframe (plot_df)
counts = plot_df.groupby(['cell_type', 'Condition'], observed=True).size().unstack(fill_value=0)
pcts = (counts / counts.sum() * 100).melt(ignore_index=False).reset_index()

# 3. Determine the sort order (Big clusters first)
cell_type_order = pcts.groupby('cell_type')['value'].sum().sort_values(ascending=False).index

# 4. Plot
plt.figure(figsize=(10, 4))
sns.barplot(
    data=pcts, 
    x='cell_type', 
    y='value', 
    hue='Condition', 
    order=cell_type_order, 
    palette={
        'Control': 'black',
        'A53T-E3': 'orange'
    },
    edgecolor='black',
    linewidth=0.5
)

plt.ylabel('Percentage of Condition (%)')
plt.xlabel('Cell Type')
plt.title('Cell Type Composition (Control vs A53T-E3)')
plt.xticks(rotation=90)
plt.legend(title='Genotype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

logging.info("Barplot generated successfully.")
plt.show()

ax = sc.pl.correlation_matrix(adata, "cell_type", figsize=(8, 8))
