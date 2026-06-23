# =============================================================================
# utils.py
# Core utility functions for snRNA-seq preprocessing and annotation.
# Imported by: 03_qc_preprocessing_annotation.py
#
# Author: Dimitrios Kyriakis
# =============================================================================

import pegasus as pg
import scanpy as sc
import anndata as ad
from anndata.tests.helpers import assert_equal
from anndata._core.sparse_dataset import SparseDataset
from anndata.experimental import read_elem, write_elem

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse
import h5py

# sys
import gc
from pathlib import Path

# pge
import sys

# pge (Pegasus Extras) must be on your PYTHONPATH before running.
# Add its directory to PYTHONPATH in your shell:
#   export PYTHONPATH="/path/to/pge:$PYTHONPATH"
# or install it if a package is available from the study authors.
import pge


def read_10x_pegasus(inpath,outpath,sample):
    h5ad_input = inpath+sample+'/outs/filtered_feature_bc_matrix.h5'
    h5ad_output = outpath+sample+'.h5ad'
    data = pg.read_input(h5ad_input)
    data.var['gene_symbol']=data.var.index
    
    data.obs['Sample'] = sample
    data.obs['Channel'] = data.obs['Sample'].copy()
    data.obs['individualID'] = data.obs['Sample'].copy()
    data.obs['cell_id'] = data.obs_names.copy()
    pge.save(data, h5ad_output)


def qc(input, output, prefix, mitocarta_file, args_mad_k=4, args_pct_mito=5, args_min_n_cells=50):
    """
    Parameters
    ----------
    mitocarta_file : str
        Path to Human.MitoCarta3.0.csv (download from https://www.broadinstitute.org/mitocarta)
    """
    data = pg.read_input(input, file_type='h5ad', genome='GRCh38', modality='rna')

    ### create channel
    data.obs['Channel'] = data.obs['Sample'].copy()
    data.obs['individualID'] = data.obs['Sample'].copy()
    data.var['gene_chrom'] = data.var['Chromosome']
    ##################
    ### QC by gene ###
    ##################

    ### identify robust genes
    pg.identify_robust_genes(data, percent_cells=0.05)

    ### remove features that are not robust (expressed at least 0.05% of cells) from downstream analysis
    data._inplace_subset_var(data.var['robust'])

    ### add ribosomal genes
    data.var['ribo'] = [x.startswith("RP") for x in data.var.gene_name]

    ### add mitochondrial genes
    data.var['mito'] = [x.startswith("MT-") for x in data.var.gene_name]

    ### add protein_coding genes
    data.var['protein_coding'] = [x == 'protein_coding' for x in data.var.gene_type]

    ### define mitocarta_genes
    mitocarta = pd.read_csv(mitocarta_file)
    data.var['mitocarta'] = [True if x in list(mitocarta.Symbol) else False for x in data.var.index]

    ### define robust_protein_coding genes (exclude ribosomal (RPL,RPS), mitochondrial, or mitocarta genes
    data.var['robust_protein_coding'] = data.var['robust'] & data.var['protein_coding']
    data.var.loc[data.var.ribo, 'robust_protein_coding'] = False
    data.var.loc[data.var.mito, 'robust_protein_coding'] = False
    data.var.loc[data.var.mitocarta, 'robust_protein_coding'] = False

    ### define robust_protein_coding_autosome genes (exclude ribosomal (RPL,RPS), mitochondrial, or mitocarta genes
    data.var['robust_protein_coding_autosome'] = data.var['robust_protein_coding'] & data.var.gene_chrom.isin(
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22'])

    ####################
    ### QC by counts ###
    ####################

    pg.qc_metrics(data, mito_prefix='MT-')

    ### nUMI and nGene QCs
    n_counts_lower, n_counts_upper = pge.qc_boundary(data.obs.n_counts, k=args_mad_k)
    print('n_UMIs lower: %d upper: %d' % (n_counts_lower, n_counts_upper))

    n_genes_lower, n_genes_upper = pge.qc_boundary(data.obs.n_genes, k=args_mad_k)
    print('n_genes lower: %d upper: %d' % (n_genes_lower, n_genes_upper))

    # log(n_counts)
    with rc_context({'figure.figsize': (4, 4)}):
        plt.figure(figsize=(4, 4))
        sns.histplot(np.log10(data.obs.n_counts))
        plt.axvline(np.log10(n_counts_lower), color='red')
        plt.axvline(np.log10(n_counts_upper), color='red')
        plt.xlabel('log10(n_counts)', fontsize=12)
        plt.savefig(prefix + "_histplot_n_counts.png")

    # log(n_genes)
    with rc_context({'figure.figsize': (4, 4)}):
        plt.figure(figsize=(4, 4))
        sns.histplot(np.log10(data.obs.n_genes))
        plt.axvline(np.log10(n_genes_lower), color='red')
        plt.axvline(np.log10(n_genes_upper), color='red')
        plt.xlabel('log10(n_genes)', fontsize=12)
        plt.savefig(prefix + "_histplot_n_genes.png")

    # scatter
    with rc_context({'figure.figsize': (4, 4)}):
        plt.figure(figsize=(4, 4))
        sns.scatterplot(x=data.obs.n_genes, y=data.obs.n_counts, alpha=0.5, s=0.1)
        plt.axhline(n_counts_lower, color='red')
        plt.axhline(n_counts_upper, color='red')
        plt.axvline(n_genes_lower, color='red')
        plt.axvline(n_genes_upper, color='red')
        plt.savefig(prefix + "_scatterplot_threshold.png")

    # percent_mito
    with rc_context({'figure.figsize': (4, 4)}):
        plt.figure(figsize=(4, 4))
        sns.histplot(data.obs.percent_mito)
        plt.axvline(args_pct_mito, color='red')
        plt.xlabel('percent_mito', fontsize=12)
        plt.savefig(prefix + "_histplot_percent_mito.png")

    ## apply QC filter
    pg.qc_metrics(data,
                  min_genes=n_genes_lower, max_genes=n_genes_upper,
                  min_umis=n_counts_lower, max_umis=n_counts_upper,
                  mito_prefix='MT-', percent_mito=args_pct_mito)

    # df = pg.get_filter_stats(data)
    # df.to_csv(prefix + "_filter_stats.csv")

    #####################
    ### QC by n_cells ###
    #####################

    n_cells_before_qc = data.obs.Channel.value_counts().rename_axis('Channel').reset_index(name='counts')
    n_cells_after_qc = data.obs[data.obs.passed_qc].Channel.value_counts().rename_axis('Channel').reset_index(
        name='counts')

    print('n_cells before QC:', np.sum(n_cells_before_qc[n_cells_before_qc.counts > 0].counts))
    print('n_cells after QC:', np.sum(n_cells_after_qc[n_cells_after_qc.counts > 0].counts))

    print('mean n_cells before QC', np.mean(n_cells_before_qc[n_cells_before_qc.counts > 0].counts))
    print('mean n_cells after QC', np.mean(n_cells_after_qc[n_cells_after_qc.counts > 0].counts))

    with rc_context({'figure.figsize': (4, 4)}):
        plt.figure(figsize=(4, 4))
        sns.histplot(np.log10(n_cells_before_qc[n_cells_before_qc.counts > 0].counts))
        plt.axvline(np.log10(args_min_n_cells), color='red')
        plt.xlabel('log10(n_cells)', fontsize=12)
        plt.savefig(prefix + "_histplot_n_cells.png")

    ### n_cells QC
    n_cells_outlier = list(n_cells_after_qc[n_cells_after_qc.counts < args_min_n_cells].Channel)
    data.obs.loc[data.obs.Channel.isin(n_cells_outlier), 'passed_qc'] = False
    print(
        'remove %i donors that have cells less than %i: %s' % (len(n_cells_outlier), args_min_n_cells, n_cells_outlier))

    ### filter cells
    pg.filter_data(data)

    ### clean unused categories
    data.obs['Channel'] = data.obs.Channel.cat.remove_unused_categories()

    ### save
    data.to_anndata().write(output)


def log1p_norm(input, output):
    adata = sc.read_h5ad(input)
    adata.layers["counts"] = adata.X

    # shifted log1p transform
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.X = sc.pp.log1p(scales_counts["X"], copy=True)

    adata.write(output)


def sig_score(input, output):
    data = pg.read_input(input)

    pg.calc_signature_score(data,
                            'cell_cycle_human')  ## 'cycle_diff', 'cycling', 'G1/S', 'G2/M' ## cell cycle gene score based on [Tirosh et al. 2015 | https://science.sciencemag.org/content/352/6282/189]
    # pg.calc_signature_score(data, 'gender_human')  # female_score, male_score
    pg.calc_signature_score(data,
                            'mitochondrial_genes_human')  # 'mito_genes' contains 13 mitocondrial genes from chrM and 'mito_ribo' contains mitocondrial ribosomal genes that are not from chrM
    pg.calc_signature_score(data, 'ribosomal_genes_human')  # ribo_genes
    pg.calc_signature_score(data, 'apoptosis_human')  # apoptosis

    pge.save(data, output)



def add_metadata(input, output, prefix, ensembl_metadata_file, mitocarta_file):
    """
    Parameters
    ----------
    ensembl_metadata_file : str
        Path to human Ensembl metadata CSV with columns 'Gene name', 'Gene stable ID',
        'Chromosome/scaffold name', 'Gene type'.
        Download from https://www.ensembl.org/biomart (GRCh38, latest release).
    mitocarta_file : str
        Path to Human.MitoCarta3.0.csv (download from https://www.broadinstitute.org/mitocarta)
    """
    data = pg.read_input(input, file_type='h5ad', genome='GRCh38', modality='rna')
    adata = data.to_anndata()
    # --- 3. Identify the common key columns based on your provided structure ---
    ensembl_metadata = pd.read_csv(human_ensembl_metadata_file)
    print(f"\nLoaded Ensembl metadata with shape: {ensembl_metadata.shape}")
    print("Ensembl metadata columns (first 5):", ensembl_metadata.columns.tolist()[:5])
    print("Ensembl metadata head:")
    print(ensembl_metadata.head())
    gene_symbol_col_in_metadata = 'Gene name'  # This is the gene symbol column
    ensembl_id_col_in_metadata = 'Gene stable ID'  # This is the Ensembl ID column
    
    # Ensure these crucial columns exist in the loaded metadata DataFrame
    if gene_symbol_col_in_metadata not in ensembl_metadata.columns:
        raise ValueError(f"Error: Gene symbol column '{gene_symbol_col_in_metadata}' not found in metadata.")
    if ensembl_id_col_in_metadata not in ensembl_metadata.columns:
        raise ValueError(f"Error: Ensembl ID column '{ensembl_id_col_in_metadata}' not found in metadata.")
    
    # --- 4. Prepare the metadata for merging ---
    ensembl_metadata_dedup = ensembl_metadata.drop_duplicates(subset=[gene_symbol_col_in_metadata], keep='first')
    ensembl_metadata_dedup = ensembl_metadata_dedup.set_index(gene_symbol_col_in_metadata)
    
    # --- 5. Align and Merge/Join the metadata to adata.var ---
    adata_var_df = adata.var.copy()
    merged_var = pd.merge(
        adata_var_df,
        ensembl_metadata_dedup,
        left_index=True,  # Use adata.var's current index (your gene symbols)
        right_index=True,  # Use the gene symbols we set as index in metadata_dedup
        how='left'  # Keep all genes from your AnnData, add matching metadata
    )
    
    # --- 6. Assign the merged DataFrame back to adata.var ---
    adata.var = merged_var
    
    print("\n--- Annotation Complete ---")
    print(f"Updated adata.var columns (first 5): {adata.var.columns.tolist()[:5]} and more...")
    print("adata.var head after annotation:")
    print(adata.var.head())
    
    # Verify some added columns
    if ensembl_id_col_in_metadata in adata.var.columns:
        print(
            f"\nNumber of genes with Ensembl IDs after merge: {adata.var[ensembl_id_col_in_metadata].count()} out of {adata.n_vars}")
        print(f"First 5 Ensembl IDs: {adata.var[ensembl_id_col_in_metadata].head().tolist()}")
    else:
        print(f"\nNote: Column '{ensembl_id_col_in_metadata}' was not found/merged into adata.var.")
    
    if 'Gene type' in adata.var.columns:  # Assuming 'Gene type' is a column you want to verify
        print(f"\nNumber of genes with Gene type: {adata.var['Gene type'].count()} out of {adata.n_vars}")
        print("Example: Check some specific genes' metadata:")
        # This will show all metadata for the first gene in your AnnData
        if not adata.var_names.empty:
            print(adata.var.loc[adata.var_names[0]])
        else:
            print("No genes left in AnnData object to display example metadata.")
    else:
        print(f"\nNote: Column 'Gene type' was not found/merged into adata.var.")
    
    # Assign standardized column names and drop original Ensembl columns
    adata.var['gene_symbol'] = adata.var.index.copy()
    adata.var['gene_name'] = adata.var.index.copy()
    # Check if 'Gene type' column exists before accessing it
    if 'Gene type' in adata.var.columns:
        adata.var['gene_type'] = adata.var['Gene type'].copy()
    else:
        print("Warning: 'Gene type' column not found in adata.var after merge. 'gene_type' will be NaN.")
        adata.var['gene_type'] = None  # Or np.nan or pd.NA
    
    # Check if 'Gene stable ID' column exists before accessing it
    if 'Gene stable ID' in adata.var.columns:
        adata.var['gene_id'] = adata.var['Gene stable ID'].copy()
    else:
        print("Warning: 'Gene stable ID' column not found in adata.var after merge. 'gene_id' will be NaN.")
        adata.var['gene_id'] = None
    
    # Check if 'Chromosome/scaffold name' column exists before accessing it
    if 'Chromosome/scaffold name' in adata.var.columns:
        adata.var['Chromosome'] = adata.var['Chromosome/scaffold name'].copy()
    else:
        print(
            "Warning: 'Chromosome/scaffold name' column not found in adata.var after merge. 'Chromosome' will be NaN.")
        adata.var['Chromosome'] = None
    
    # Drop original Ensembl columns if they exist
    cols_to_drop = []
    if 'Gene type' in adata.var.columns:
        cols_to_drop.append('Gene type')
    if 'Gene stable ID' in adata.var.columns:
        cols_to_drop.append('Gene stable ID')
    if 'Chromosome/scaffold name' in adata.var.columns:
        cols_to_drop.append('Chromosome/scaffold name')
    
    if cols_to_drop:
        adata.var.drop(columns=cols_to_drop, inplace=True)
    data.var = adata.var
    mitocarta = pd.read_csv(mitocarta_file)
    data.var['mitocarta'] = [True if x in list(mitocarta.Symbol) else False for x in data.var.index]
    # Compute total counts per cell
    data.obs["n_counts"] = np.array(data.X.sum(axis=1)).flatten()
    
    # Compute number of detected genes per cell (nonzero features)
    data.obs["n_features"] = np.array((data.X > 0).sum(axis=1)).flatten()

    pge.save(data, output)

#!/usr/bin/env python
"""
Global clustering and UMAP embedding for large snRNA-seq datasets.

This script demonstrates how to:
  1) Subset an AnnData .h5ad file on disk to protein-coding autosomal genes,
  2) Select highly variable genes (HVGs),
  3) Run PCA + Harmony batch correction,
  4) Build a kNN graph and compute UMAP,
  5) Run Leiden clustering at two resolutions (0.10 and 0.20),
  6) Save the resulting AnnData object.

Dependencies:
  - pegasus            (https://github.com/lilab-bcb/pegasus)
  - scanpy             (https://scanpy.readthedocs.io)
  - pge (Pegasus extras; optional, or replace with your own utilities)

This script assumes:
  - Input is a large .h5ad file with raw counts in X,
  - .var contains: gene_type, gene_chrom (chromosome annotation),
  - .obs contains: a batch key (e.g. "Brain_bank" or "Source").

Author: Tereza Clarence, Donghoon Lee

### Global clustering & UMAP

Example usage:

```bash
python scripts/global_clustering_umap.py \
  --orig-h5ad data/AMP_PD_freeze2_all.h5ad \
  --subset-h5ad data/AMP_PD_freeze2_autosome_pc.h5ad \
  --output-h5ad results/AMP_PD_freeze2_global_umap.h5ad \
  --batch-key Brain_bank
"""

import argparse
from pathlib import Path

import pegasus as pg
import scanpy as sc
import h5py

import numpy as np

# If pge is part of your repo, import it; otherwise replace with your own utilities
import pge  # expected to provide read_everything_but_X() and ondisk_subset(), scanpy_hvf_h5ad()


def subset_on_disk(
    orig_h5ad: str,
    new_h5ad: str,
    chunk_size: int = 500_000,
):
    """
    Create an on-disk subset of an AnnData file, keeping only protein-coding autosomal genes.
    """
    print(f"[INFO] Reading metadata from {orig_h5ad}")
    adata = pge.read_everything_but_X(orig_h5ad)

    # Keep all cells
    subset_obs = (adata.obs_names != None).tolist()

    # Keep autosomal protein-coding genes (exclude MT, X, Y)
    subset_var = (
        (adata.var["gene_type"] == "protein_coding")
        & (~adata.var["gene_chrom"].isin(["MT", "X", "Y"]))
    ).tolist()

    print("[INFO] Subsetting on disk to protein-coding autosomal genes...")
    pge.ondisk_subset(
        orig_h5ad=orig_h5ad,
        new_h5ad=new_h5ad,
        subset_obs=subset_obs,
        subset_var=subset_var,
        chunk_size=chunk_size,
        raw=False,
    )

    # Optional: quick structural check
    with h5py.File(new_h5ad, "r") as f:
        print("[INFO] New X dataset structure:")
        f["X"].visititems(print)

    print(f"[DONE] Wrote subsetted file to {new_h5ad}")


def run_global_clustering_umap(
    input_h5ad: str,
    output_h5ad: str,
    batch_key: str = "Brain_bank",
    n_top_genes: int = 6000,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
    n_pcs: int = 30,
    k_graph: int = 100,
    n_neighbors_umap: int = 100,
):
    """
    Perform global clustering and UMAP on a (potentially large) snRNA-seq dataset.
    """

    # -------------------------------------------------------------------------
    # 1. Select HVGs using Scanpy-style HVG function
    # -------------------------------------------------------------------------
    print(f"[INFO] Selecting HVGs from {input_h5ad}")
    hvg = pge.scanpy_hvf_h5ad(
        h5ad_file=input_h5ad,
        flavor="cell_ranger",
        batch_key=batch_key,
        n_top_genes=n_top_genes,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        protein_coding=True,
        autosome=True,
    )
    print(f"[INFO] Selected {len(hvg)} HVGs")

    # Optional: save HVG list
    hvg_txt = Path(output_h5ad).with_suffix(".hvg.txt")
    with open(hvg_txt, "w") as fp:
        for g in hvg:
            fp.write(f"{g}\n")
    print(f"[INFO] Saved HVG list to {hvg_txt}")

    # Re-read to ensure consistent ordering
    with open(hvg_txt) as f:
        hvg = [line.rstrip() for line in f]

    # -------------------------------------------------------------------------
    # 2. Load data via Pegasus and mark HVGs
    # -------------------------------------------------------------------------
    print(f"[INFO] Loading data from {input_h5ad}")
    data = pg.read_input(input_h5ad, genome="GRCh38", modality="rna")
    print(data)

    data.var["highly_variable_features"] = False
    data.var.loc[data.var.index.isin(hvg), "highly_variable_features"] = True

    print("[INFO] HVG flag summary:")
    print(data.var["highly_variable_features"].value_counts())
    print("[INFO] HVGs by chromosome:")
    print(data.var[data.var["highly_variable_features"] == True]["gene_chrom"].value_counts())

    # -------------------------------------------------------------------------
    # 3. PCA, regression, Harmony integration
    # -------------------------------------------------------------------------
    print("[INFO] Running PCA")
    pg.pca(data, n_components=n_pcs)
    pg.elbowplot(data)
    npc = min(data.uns["pca_ncomps"], n_pcs)
    print(f"[INFO] Using {npc} PCs for downstream steps")

    print("[INFO] Regressing out technical covariates and running Harmony")
    pg.regress_out(data, attrs=["n_counts", "percent_mito", "cycle_diff"])
    pg.run_harmony(
        data,
        batch=batch_key,
        rep="pca_regressed",
        max_iter_harmony=20,
        n_comps=npc,
    )

    print("[INFO] Building kNN graph")
    pg.neighbors(
        data,
        rep="pca_regressed_harmony",
        use_cache=False,
        dist="l2",
        K=k_graph,
        n_comps=npc,
    )

    # -------------------------------------------------------------------------
    # 4. UMAP and Leiden clustering
    # -------------------------------------------------------------------------
    print("[INFO] Computing UMAP")
    pg.umap(
        data,
        rep="pca_regressed_harmony",
        n_neighbors=n_neighbors_umap,
        rep_ncomps=npc,
    )

    adata = data.to_anndata()

    print("[INFO] Running Leiden clustering (res 0.10, 0.20)")
    sc.pp.neighbors(adata, use_rep="X_pca_regressed_harmony")
    sc.tl.leiden(adata, key_added="leiden_res0_10", resolution=0.10)
    sc.tl.leiden(adata, key_added="leiden_res0_20", resolution=0.20)

    # Save final AnnData
    adata.write_h5ad(output_h5ad)
    print(f"[DONE] Saved clustered AnnData to {output_h5ad}")

    # Optional: UMAP figure coloured by Leiden clusters
    sc.pl.umap(
        adata,
        color=["leiden_res0_20", "leiden_res0_10"],
        legend_loc="on data",
        frameon=False,
        legend_fontsize=5,
        legend_fontoutline=1,
        title=["Leiden res 0.20", "Leiden res 0.10"],
        size=1,
        wspace=0,
        ncols=2,
        save="_global_clustering.png",
    )
    print("[DONE] UMAP plot saved (scanpy default location)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Global clustering + UMAP for large snRNA-seq AnnData files."
    )
    parser.add_argument(
        "--orig-h5ad",
        type=str,
        required=True,
        help="Path to the original (large) AnnData .h5ad file.",
    )
    parser.add_argument(
        "--subset-h5ad",
        type=str,
        required=True,
        help="Path for the on-disk subset (protein-coding autosomal genes).",
    )
    parser.add_argument(
        "--output-h5ad",
        type=str,
        required=True,
        help="Path for the final clustered .h5ad file (with UMAP + Leiden).",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default="Brain_bank",
        help="Batch key used for Harmony integration (e.g. brain bank / source).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500_000,
        help="Chunk size for on-disk subsetting.",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=6000,
        help="Number of HVGs to select globally.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of principal components to retain.",
    )
    parser.add_argument(
        "--k-graph",
        type=int,
        default=100,
        help="Number of neighbors (K) for kNN graph construction.",
    )
    parser.add_argument(
        "--n-neighbors-umap",
        type=int,
        default=100,
        help="Number of neighbors for UMAP.",
    )
    return parser.parse_args()


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
    
    # Optional: Rotate row labels if they are long
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0) 
    plt.show()

    if legend_on:
        sc.pl.umap(adata,color=[cluster,agr_name+'_corr_annot_pnas'],legend_loc='on data')
    else:
        sc.pl.umap(adata,color=[cluster,agr_name+'_corr_annot_pnas'])
    return(corr_mat_pge)

