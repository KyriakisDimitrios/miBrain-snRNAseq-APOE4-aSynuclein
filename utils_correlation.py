# =============================================================================
# utils_correlation.py
# Cross-dataset pseudobulk correlation utilities.
# Imported by: 04_annotation_concordance_AMPPD.py
#
# Author: Dimitrios Kyriakis
# =============================================================================

import os
# sc
import pegasus as pg
import scanpy as sc
import anndata as ad

# data
import numpy as np
import pandas as pd

# plot
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns
from adjustText import adjust_text

sns.set_style("whitegrid", {'axes.grid': False})
sc.set_figure_params(scanpy=True, dpi=100, dpi_save=300, fontsize=12, color_map='YlOrRd')  # 'viridis_r'
sc.settings.verbosity = 1
sc.logging.print_header()

# sys
import gc
from pathlib import Path
import time
import logging

logger = logging.getLogger("pegasus")

import logging
import sys
import random

import numpy as np
import scanpy as sc
import pegasus as pg

# plotting
import matplotlib.pyplot as plt
import seaborn as sns

from concurrent.futures import ProcessPoolExecutor
import warnings
import pge

def get_nominal_degs_fast(h5ad_file, cluster_label, p_thresh=0.05, mat_key='X',
                          subset=True, max_cells_per_group=500, seed=42):
    """
    Subsamples per cluster (stratified) and extracts nominal DEGs with optimized memory handling.
    """
    warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
    
    # 1. Load and Filter
    adata = sc.read_h5ad(h5ad_file)
    counts = adata.obs[cluster_label].value_counts()
    keep_clusters = counts[counts >= 2].index.tolist()
    adata = adata[adata.obs[cluster_label].isin(keep_clusters)].copy()

    # 2. Stratified Subsampling (Corrected Logic)
    if subset:
        logging.info(f"Stratified subsampling: max {max_cells_per_group} cells per group in {cluster_label}")
        np.random.seed(seed)
        sampled_indices = []
        
        for cluster in keep_clusters:
            cluster_indices = adata.obs_names[adata.obs[cluster_label] == cluster]
            if len(cluster_indices) > max_cells_per_group:
                choice = np.random.choice(cluster_indices, max_cells_per_group, replace=False)
                sampled_indices.extend(choice)
            else:
                sampled_indices.extend(cluster_indices)
        
        adata = adata[sampled_indices].copy()

    # 3. Prepare matrix for DEG computation
    if mat_key == 'counts':
        adata.X = adata.layers[mat_key].copy()
    else:
        adata.layers['counts'] = adata.X.copy()
        
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 4. Compute DEGs (Wilcoxon)
    sc.tl.rank_genes_groups(adata, cluster_label, method='wilcoxon',use_raw=False)
    
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_sig_genes = set()
    for group in groups:
        mask = result['pvals'][group] < p_thresh
        all_sig_genes.update(result['names'][group][mask])
        
    logging.info(f"Extracted {len(all_sig_genes)} nominal DEGs from {len(groups)} clusters.")
    return list(all_sig_genes)
    
def get_fdr_degs_fast(h5ad_file, cluster_label, p_thresh=0.05, mat_key='X',
                          subset=True, max_cells_per_group=500, seed=42):
    """
    Subsamples per cluster (stratified) and extracts FDR-significant DEGs with optimized memory handling.
    """
    warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
    
    # 1. Load and Filter
    adata = sc.read_h5ad(h5ad_file)
    counts = adata.obs[cluster_label].value_counts()
    keep_clusters = counts[counts >= 2].index.tolist()
    adata = adata[adata.obs[cluster_label].isin(keep_clusters)].copy()

    # 2. Stratified Subsampling
    if subset:
        logging.info(f"Stratified subsampling: max {max_cells_per_group} cells per group in {cluster_label}")
        np.random.seed(seed)
        sampled_indices = []
        
        for cluster in keep_clusters:
            cluster_indices = adata.obs_names[adata.obs[cluster_label] == cluster]
            if len(cluster_indices) > max_cells_per_group:
                choice = np.random.choice(cluster_indices, max_cells_per_group, replace=False)
                sampled_indices.extend(choice)
            else:
                sampled_indices.extend(cluster_indices)
        
        adata = adata[sampled_indices].copy()

    # 3. Prepare matrix for DEG computation
    if mat_key == 'counts':
        adata.X = adata.layers[mat_key].copy()
    else:
        adata.layers['counts'] = adata.X.copy()
        
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 4. Compute DEGs (Wilcoxon)
    sc.tl.rank_genes_groups(adata, cluster_label, method='wilcoxon',use_raw=False)
    
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_sig_genes = set()
    for group in groups:
        mask = result['pvals_adj'][group] < p_thresh  # FDR correction applied here
        all_sig_genes.update(result['names'][group][mask])
        
    logging.info(f"Extracted {len(all_sig_genes)} FDR-significant DEGs from {len(groups)} clusters.")
    return list(all_sig_genes)


def preprocess_for_correlation(df, scale="standard", min_variance_percentile=50):
    """
    Prepares pseudobulk matrix for correlation:
    - Removes zero-variance genes
    - Filters out low-variance (housekeeping) genes
    - scale="robust": Uses Median/MAD to protect against outlier-driven skew
    - scale="standard": Uses Mean/Std
    """
    # 1. Calculate Standard Deviation once to save memory/time
    gene_std = df.std(axis=0)
    
    # 2. Remove absolute flatliners
    df = df.loc[:, gene_std > 0]
    if df.empty:
        logging.warning("Matrix empty after zero-variance removal.")
        return df

    # 3. Variance Percentile Filter (Keep only the 'strong heartbeats')
    # 
    gene_std = gene_std[gene_std > 0]
    var_threshold = np.percentile(gene_std, min_variance_percentile)
    df = df.loc[:, gene_std >= var_threshold]
    
    logging.info(f"Genes after variance filtering (>{min_variance_percentile}th pct): {df.shape[1]}")

    # 4. Scaling Logic
    if scale == "robust":
        # Robust Z-score: (x - median) / (MAD * 1.4826)
        # 
        median = df.median(axis=0)
        mad = (df - median).abs().median(axis=0)
        
        # Avoid division by zero for genes with identical values across clusters
        mad = mad.replace(0, np.nan)
        
        # 1.4826 makes MAD consistent with the standard deviation of a normal distribution
        df = (df - median) / (mad * 1.4826)
        
        # Clean up any infinites caused by extremely low MAD
        pd.set_option('mode.use_inf_as_na', True)
        df.dropna(axis=1, how="any", inplace=True)
        pd.reset_option('mode.use_inf_as_na')
        logging.info("Applied Robust Scaling (Median/MAD).")
    elif scale == "standard":
        df = (df - df.mean(axis=0)) / df.std(axis=0)
        df.dropna(axis=1, how="any", inplace=True)
        logging.info("Applied Standard Scaling (Mean/Std).")
    else:
        logging.info("No Scaling,")
    return df.sort_index(axis=1)
    

def plot_correlation(
    final_corr_df,
    method='spearman',
    prefix='prefix',
    show_plot=True,
    save_plot=False,
    figure_path='./',
    annot=True,
    fontsize_annot=9,
    fontsize_ticks=10,
    fontsize_title=13,
    fontsize_cbar=10,
    fig_width=None,
    fig_height=None,
    x_title=None,
    y_title=None,
    show_title=True   # set False to hide main title
):
    plot_df = final_corr_df.dropna(how='any', axis=0)
    plot_df = plot_df.dropna(how='any', axis=1)

    if plot_df.empty:
        logging.warning(f"Empty matrix for {method}. Skipping.")
        return

    vals      = plot_df.values.flatten()
    vals      = vals[~np.isnan(vals)]
    vmin_data = np.percentile(vals, 2)
    vmax_data = np.percentile(vals, 98)

    has_negatives = vmin_data < -0.05
    has_positives = vmax_data >  0.05
    is_diverging  = has_negatives and has_positives

    if is_diverging:
        max_abs = max(abs(vmin_data), abs(vmax_data))
        vmin    = -max_abs
        vmax    =  max_abs
        center  = 0
        cmap    = "vlag"
        logging.info(f"{method}: diverging | vmin={vmin:.3f} vmax={vmax:.3f}")
    else:
        vmin   = vmin_data
        vmax   = vmax_data
        center = None
        cmap   = "Reds" if has_positives else "YlGnBu_r"
        logging.info(f"{method}: sequential ({cmap}) | vmin={vmin:.3f} vmax={vmax:.3f}")

    num_rows, num_cols = plot_df.shape
    cell_size          = max(0.3, min(0.6, 8 / max(num_rows, num_cols)))
    dynamic_height     = max(4, num_rows * cell_size + 2)
    dynamic_width      = max(6, num_cols * cell_size + 3)

    figsize = (
        fig_width  if fig_width  is not None else dynamic_width,
        fig_height if fig_height is not None else dynamic_height
    )
    logging.info(f"Figure size: {figsize[0]:.1f} x {figsize[1]:.1f} inches "
                 f"({'manual' if fig_width or fig_height else 'dynamic'})")

    g = sns.clustermap(
        plot_df,
        metric    = "correlation",
        method    = "average",
        cmap      = cmap,
        center    = center,
        vmin      = vmin,
        vmax      = vmax,
        figsize   = figsize,
        annot     = annot,
        fmt       = ".2f",
        annot_kws = {"size": fontsize_annot}
    )

    g.ax_cbar.set_ylabel(f'{method} r', fontsize=fontsize_cbar)
    g.ax_cbar.yaxis.set_label_position('left')

    if show_title:
        scale_note = "diverging | white=0" if is_diverging else f"sequential | range [{vmin:.2f}–{vmax:.2f}]"
        plt.suptitle(
            f"{prefix}\n{method} correlation",
            fontsize=fontsize_title, y=1.02
        )

    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0,  fontsize=fontsize_ticks)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, fontsize=fontsize_ticks, ha='right')

    if x_title is not None:
        g.ax_heatmap.set_xlabel(x_title, fontsize=fontsize_ticks + 1, fontweight='bold', labelpad=10)
    else:
        g.ax_heatmap.set_xlabel('')

    if y_title is not None:
        g.ax_heatmap.set_ylabel(y_title, fontsize=fontsize_ticks + 1, fontweight='bold', labelpad=10)
    else:
        g.ax_heatmap.set_ylabel('')

    if save_plot:
        os.makedirs(figure_path, exist_ok=True)
        save_name = os.path.join(figure_path, f"{prefix}_{method}_correlation.png")
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
        save_name = os.path.join(figure_path, f"{prefix}_{method}_correlation.pdf")
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
        logging.info(f"Saved: {save_name}")

    if show_plot:
        plt.show()
    else:
        plt.close()

        
def correlation_qc_report(
    spearman_df,
    pearson_df,
    degs_h5ad1,
    degs_h5ad2,
    common_sig_genes,
    prefix='prefix',
    save_report=False,
    figure_path='./'
):
    print("\n" + "=" * 60)
    print(f"  CORRELATION QC REPORT — {prefix}")
    print("=" * 60)

    print(f"\n[DEG Summary]")
    print(f"  h5ad1 DEGs (FDR<0.05):   {len(degs_h5ad1):>8}")
    print(f"  h5ad2 DEGs (FDR<0.05):   {len(degs_h5ad2):>8}")
    print(f"  Final feature set:        {len(common_sig_genes):>8}")

    for method, df in [('Spearman', spearman_df), ('Pearson', pearson_df)]:
        vals = df.values.flatten()
        vals = vals[~np.isnan(vals)]
        print(f"\n[{method} Distribution]")
        print(f"  Shape:      {df.shape}")
        print(f"  Mean:       {vals.mean():.3f}")
        print(f"  Std:        {vals.std():.3f}")
        print(f"  Min:        {vals.min():.3f}")
        print(f"  Max:        {vals.max():.3f}")
        print(f"  > 0.5:      {(vals > 0.5).sum():>6} pairs  ({100*(vals>0.5).mean():.1f}%)")
        print(f"  > 0.7:      {(vals > 0.7).sum():>6} pairs  ({100*(vals>0.7).mean():.1f}%)")
        print(f"  < 0:        {(vals < 0.0).sum():>6} pairs  ({100*(vals<0.0).mean():.1f}%)")

    print(f"\n[Best Match per Cluster — Spearman]")
    best_matches = spearman_df.idxmax(axis=1)
    best_scores  = spearman_df.max(axis=1).round(3)
    worst_scores = spearman_df.min(axis=1).round(3)
    match_df = pd.DataFrame({
        'best_match' : best_matches,
        'max_r'      : best_scores,
        'min_r'      : worst_scores,
        'delta'      : (best_scores - worst_scores).round(3)
    }).sort_values('max_r', ascending=False)
    print(match_df.to_string())
    print("\n  [delta = max_r - min_r: higher = more discriminating]")

    # Distribution histograms
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    for ax, (method, df) in zip(axes, [('Spearman', spearman_df), ('Pearson', pearson_df)]):
        vals = df.values.flatten()
        vals = vals[~np.isnan(vals)]
        ax.hist(vals, bins=30, color='steelblue', edgecolor='white', alpha=0.85)
        ax.axvline(0,             color='black',  linestyle='-',  lw=1,
                   label='r=0')
        ax.axvline(vals.mean(),   color='orange', linestyle='--', lw=1.5,
                   label=f'mean={vals.mean():.2f}')
        ax.axvline(np.median(vals), color='red',  linestyle=':',  lw=1.5,
                   label=f'median={np.median(vals):.2f}')
        ax.set_title(f'{method} Distribution', fontsize=12)
        ax.set_xlabel('Correlation coefficient')
        ax.set_ylabel('Count')
        ax.legend(fontsize=9)

    plt.suptitle(f'{prefix} — Correlation QC', fontsize=13)
    plt.tight_layout()

    if save_report:
        os.makedirs(figure_path, exist_ok=True)
        save_name = os.path.join(figure_path, f"{prefix}_QC_report.png")
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
        logging.info(f"QC report saved: {save_name}")

    plt.show()
    print("=" * 60 + "\n")


def correlation_between_h5ads(
    h5ad1, h5ad2,
    group_by1='subclass', group_by2='subclass',
    mat_key1='counts', mat_key2='counts',
    degs_calc='first',
    nominal = True,
    scale_before_spearman = False,
    scale = 'robust',  # 'standard' or False
    subset1=True, subset2=True,
    max_cells_per_group1=500,
    max_cells_per_group2=500,
    min_variance_percentile=50,
    deg_intersection=False,
    flip=False,
    prefix='prefix',
    show_plots=True,
    save_plots=False,
    figure_path='./figures/'
):
    """
    Cross-dataset cluster correlation via pseudobulk profiles.

    Statistical design:
    ┌──────────────────────────────────────────────────────────┐
    │ DEGs:        subsampled cells, Wilcoxon + BH FDR         │
    │              → feature selection only                    │
    │ Correlation: pseudobulk cluster means (all cells)        │
    │   Spearman:  rank-based, no scaling needed               │
    │   Pearson:   robust Z-score (median/MAD) for small N     │
    │ Variance filter: removes housekeeping genes              │
    │ Colorbar:    symmetric around 0, data-driven max         │
    │ NaN:         dropped, never filled                       │
    └──────────────────────────────────────────────────────────┘
    """

    # 1. Pseudobulk
    logging.info("Aggregating pseudobulk matrices...")
    agr_mat_h5ad1 = pge.pb_agg_by_cluster(h5ad1, group_by1, log1p=True, mat_key=mat_key1)
    agr_mat_h5ad2 = pge.pb_agg_by_cluster(h5ad2, group_by2, log1p=True, mat_key=mat_key2)

    genes_h5ad1_set = set(agr_mat_h5ad1.columns)
    genes_h5ad2_set = set(agr_mat_h5ad2.columns)
    common_features = genes_h5ad1_set.intersection(genes_h5ad2_set)
    logging.info(f"Shared genes: {len(common_features)}")

    # 2. DEGs
    degs_h5ad1, degs_h5ad2 = set(), set()

    if degs_calc in ['first', 'both']:
        
        if nominal:
            logging.info(f"Nominal DEGs h5ad1 (subset={subset1}, max={max_cells_per_group1})...")
            degs_h5ad1 = set(get_nominal_degs_fast(
                h5ad1, group_by1, mat_key=mat_key1,
                subset=subset1, max_cells_per_group=max_cells_per_group1
            ))
        else:
           logging.info(f"FDR DEGs h5ad1 (subset={subset1}, max={max_cells_per_group1})...")
           degs_h5ad1 = set(get_fdr_degs_fast(
                h5ad1, group_by1, mat_key=mat_key1,
                subset=subset1, max_cells_per_group=max_cells_per_group1
            )) 
        
        logging.info(f"h5ad1 DEGs: {len(degs_h5ad1)}")

    if degs_calc in ['second', 'both']:
        if nominal:
            logging.info(f"Nominal DEGs h5ad2 (subset={subset2}, max={max_cells_per_group2})...")
            degs_h5ad2 = set(get_nominal_degs_fast(
                h5ad2, group_by2, mat_key=mat_key2,
                subset=subset2, max_cells_per_group=max_cells_per_group2
            ))
        else:
            logging.info(f"FDR DEGs h5ad2 (subset={subset2}, max={max_cells_per_group2})...")
            degs_h5ad2 = set(get_fdr_degs_fast(
                h5ad2, group_by2, mat_key=mat_key2,
                subset=subset2, max_cells_per_group=max_cells_per_group2
            ))
        logging.info(f"h5ad2 DEGs: {len(degs_h5ad2)}")

    # 3. Feature selection
    if degs_calc == 'first':
        common_sig_genes = sorted(common_features.intersection(degs_h5ad1))
    elif degs_calc == 'second':
        common_sig_genes = sorted(common_features.intersection(degs_h5ad2))
    elif degs_calc == 'both':
        combined_degs = (
            degs_h5ad1.intersection(degs_h5ad2) if deg_intersection
            else degs_h5ad1.union(degs_h5ad2)
        )
        logging.info(f"DEG mode: {'intersection' if deg_intersection else 'union'}")
        common_sig_genes = sorted(common_features.intersection(combined_degs))

    logging.info(f"Final feature set: {len(common_sig_genes)} genes.")

    if len(common_sig_genes) == 0:
        logging.error("No genes remain. Relax p_thresh or check group_by columns.")
        return None, None

    # 4. Preprocessing
    # Spearman — no scaling (rank-based)
    if scale_before_spearman:
        mat1_unscaled = preprocess_for_correlation(
            agr_mat_h5ad1[common_sig_genes],
            scale=scale,
            min_variance_percentile=min_variance_percentile
        )
        mat2_unscaled = preprocess_for_correlation(
            agr_mat_h5ad2[common_sig_genes],
            scale=scale,
            min_variance_percentile=min_variance_percentile
        )
    else:
        mat1_unscaled = preprocess_for_correlation(
            agr_mat_h5ad1[common_sig_genes],
            scale=False,
            min_variance_percentile=min_variance_percentile
        )
        mat2_unscaled = preprocess_for_correlation(
            agr_mat_h5ad2[common_sig_genes],
            scale=False,
            min_variance_percentile=min_variance_percentile
        )

    # Pearson — robust Z-score (median/MAD)
    mat1_scaled = preprocess_for_correlation(
        agr_mat_h5ad1[common_sig_genes],
        scale=scale,
        min_variance_percentile=min_variance_percentile
    )
    mat2_scaled = preprocess_for_correlation(
        agr_mat_h5ad2[common_sig_genes],
        scale=scale,
        min_variance_percentile=min_variance_percentile
    )

    # Align genes after filtering
    spearman_genes = sorted(set(mat1_unscaled.columns).intersection(mat2_unscaled.columns))
    pearson_genes  = sorted(set(mat1_scaled.columns).intersection(mat2_scaled.columns))

    logging.info(f"Spearman genes: {len(spearman_genes)} | Pearson genes: {len(pearson_genes)}")

    # 5. Spearman
    spearman_df = pge.spearman_corr(
        mat1_unscaled[spearman_genes],
        mat2_unscaled[spearman_genes]
    )

    # 6. Pearson on robust Z-scored pseudobulk
    pearson_df = pge.pearson_corr(
        mat1_scaled[pearson_genes],
        mat2_scaled[pearson_genes]
    )

    # 7. Clean
    for df in [spearman_df, pearson_df]:
        df.dropna(how='all', axis=0, inplace=True)
        df.dropna(how='all', axis=1, inplace=True)

    if flip:
        spearman_df = spearman_df.T
        pearson_df  = pearson_df.T

    # 8. Plot + QC
    final_prefix = (
        f"{prefix}_DEGs_{degs_calc}_"
        f"{'intersect' if deg_intersection else 'union'}_"
        f"var{min_variance_percentile}"
    )
    plot_opts = {
        'prefix'     : final_prefix,
        'show_plot'  : show_plots,
        'save_plot'  : save_plots,
        'figure_path': figure_path
    }

    plot_correlation(spearman_df, method='spearman', **plot_opts)
    plot_correlation(pearson_df,  method='pearson',  **plot_opts)

    correlation_qc_report(
        spearman_df, pearson_df,
        degs_h5ad1, degs_h5ad2,
        common_sig_genes,
        prefix=final_prefix,
        save_report=save_plots,
        figure_path=figure_path
    )

    return spearman_df, pearson_df