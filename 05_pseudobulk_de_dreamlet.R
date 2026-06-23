# =============================================================================
# 05_pseudobulk_de_dreamlet.R
# Pseudobulk differential expression (dreamlet 1.9) and pathway enrichment
# (fgsea) across cell types in A53T miBRAIN snRNA-seq data.
#
# Author: Dimitrios Kyriakis
#
# Usage:
#   Rscript 05_pseudobulk_de_dreamlet.R \
#       --workdir    /path/to/working/directory \
#       --prefix     miBRAIN_ \
#       --kegg-gmt   /path/to/c2.cp.kegg.v7.5.1.symbols.gmt \
#       --gobp-gmt   /path/to/c5.go.bp.v7.5.1.symbols.gmt \
#       --gomf-gmt   /path/to/c5.go.mf.v7.5.1.symbols.gmt \
#       --outdir     /path/to/output
#
# Arguments:
#   --workdir    Directory containing input CSV files (counts_matrix.csv,
#                cell_metadata.csv, umap_coordinates.csv) and for outputs
#   --prefix     Filename prefix for input CSVs [default: miBRAIN_]
#   --kegg-gmt   KEGG GMT file (MSigDB v7.5.1)
#                Download: https://www.gsea-msigdb.org/gsea/downloads.jsp
#   --gobp-gmt   GO Biological Process GMT file (MSigDB v7.5.1)
#   --gomf-gmt   GO Molecular Function GMT file (MSigDB v7.5.1)
#   --outdir     Output directory for results and figures [default: --workdir]
# =============================================================================

library(optparse)

option_list <- list(
  make_option("--workdir",  type="character", default=NULL,
              help="Working directory with input CSVs"),
  make_option("--prefix",   type="character", default="miBRAIN_",
              help="Prefix for input CSV files [default: miBRAIN_]"),
  make_option("--kegg-gmt", type="character", default=NULL,
              help="Path to KEGG GMT file (MSigDB)"),
  make_option("--gobp-gmt", type="character", default=NULL,
              help="Path to GO:BP GMT file (MSigDB)"),
  make_option("--gomf-gmt", type="character", default=NULL,
              help="Path to GO:MF GMT file (MSigDB)"),
  make_option("--outdir",   type="character", default=NULL,
              help="Output directory [default: same as --workdir]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$workdir)) stop("--workdir is required")
setwd(opt$workdir)
if (is.null(opt$outdir)) opt$outdir <- opt$workdir
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

INPUT_PREFIX <- opt$prefix
KEGG_GMT     <- opt[["kegg-gmt"]]
GOBP_GMT     <- opt[["gobp-gmt"]]
GOMF_GMT     <- opt[["gomf-gmt"]]
OUTDIR       <- opt$outdir

library(zellkonverter)
library(SingleCellExperiment)
library(dreamlet)
library(ggplot2)
library(tidyverse)
library(aplot)
library(ggtree)
library(scattermore)
library(zenith)
library(crumblr)
library(GSEABase)
library(qvalue)
library(BiocParallel)
library(cowplot)
library(presto)
library(DelayedArray)

# --- Print Versions for Reproducibility ---

cat("Session package versions:\n")
cat(paste0("dreamlet v", packageVersion("dreamlet")), "\n")
cat(paste0("crumblr v", packageVersion("crumblr")), "\n")
cat(paste0("variancePartition v", packageVersion("variancePartition")), "\n")
cat(paste0("zenith v", packageVersion("zenith")), "\n")
cat(paste0("zellkonverter v", packageVersion("zellkonverter")), "\n")
cat(paste0("BiocManager v", BiocManager::version()), "\n")

print(getwd())


# --- Step 1: Load H5AD and Metadata ---

library(Seurat)
options(Seurat.object.assay.version = 'v3')

# Working directory set via --workdir argument above
getwd()

sessionInfo()

# # Read csv files and crete Seurat Object

rerun = TRUE
if (rerun){
    # Load the counts matrix and cell metadata from CSV files
    # prefix = '2026_03_02_Scanpy_pct020_annotated_'
    prefix = INPUT_PREFIX
    counts_matrix <- read.csv(paste0(prefix , "counts_matrix.csv"), row.names = 1)
    cell_metadata <- read.csv(paste0(prefix , "cell_metadata.csv"), row.names = 1)
    umap_coords <- read.csv(paste0(prefix , "umap_coordinates.csv"), row.names = 1)
    
    # Create the Seurat object
    seurat_obj <- CreateSeuratObject(counts = t(counts_matrix), meta.data = cell_metadata)
    
    # Optionally, if you have a specific assay (e.g., RNA), set it
    DefaultAssay(seurat_obj) <- "RNA"
    
    # View the Seurat object to ensure it's created correctly
    seurat_obj
    
    umap_coords_set = umap_coords[,c('UMAP_1',	'UMAP_2')]
    seurat_obj[["umap"]] <- CreateDimReducObject(
      embeddings = as.matrix(umap_coords_set), # Must be a matrix
      key = "UMAP_", # Key must end with an underscore, e.g., UMAP_1, UMAP_2
      assay = DefaultAssay(seurat_obj) # Assign to the correct assay
    )
        
    rm(cell_metadata)
    rm(counts_matrix)
    rm(umap_coords)
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    saveRDS(seurat_obj,'miBRAIN_Seurat.rds')
}

# # Load Seurat Object

library(Seurat)
seurat_obj <- readRDS('miBRAIN_Seurat.rds')

seurat_obj$class <- seurat_obj$cell_type

# Extract as character to prevent factor-level mismatch errors
conditions <- as.character(seurat_obj$Condition)

# 1. Exact match replacement for 'wt' -> 'Cntrl' (logical indexing is highly optimized)
conditions[conditions == "wt"] <- "Cntrl"

# 2. String replacement for hyphens -> underscores
conditions <- gsub("-", "_", conditions)

# Reassign back to the Seurat object
seurat_obj$dx <- conditions

as.data.frame(table(seurat_obj$cell_type,seurat_obj$dx))

write_csv(as.data.frame(table(seurat_obj$cell_type,seurat_obj$dx)),'Number_of_cells_per_dx_ct')

write_csv(as.data.frame(table(seurat_obj$cell_type,seurat_obj$Sample)),'Number_of_cells_per_sample_ct')

write.csv(seurat_obj@meta.data,'Metadata.csv')

library(dittoSeq)
dittoBarPlot(seurat_obj, "class", group.by = "Sample")
dittoBarPlot(seurat_obj, "class", group.by = "Condition")

DimPlot(seurat_obj, reduction = "umap",group.by = 'leiden_res0_20',label=T,raster=T)

DimPlot(seurat_obj, reduction = "umap",group.by = 'class',label=T,raster=T)+NoLegend()

library(Seurat)
library(presto)   # auto-detected by Seurat — makes Wilcoxon 10-100x faster
library(dplyr)

Idents(seurat_obj) <- "class"

# ── Stage 1: FindAllMarkers with specificity filters ─────────────────────────
# min.diff.pct=0.2 enforces cluster specificity AND speeds up the test
# by skipping genes that don't meet the spread criterion before Wilcoxon runs
cat(sprintf("[%s] INFO | Running FindAllMarkers...\n", Sys.time()), file = stderr())

fast_markers <- FindAllMarkers(
  seurat_obj,
  only.pos        = TRUE,
  test.use        = "wilcox",
  logfc.threshold = 0.5,
  min.pct         = 0.2,
  min.diff.pct    = 0.2,
  return.thresh   = 0.05   # nominal p-value pre-filter to reduce output size
)

# Save raw output — always keep this, useful for re-filtering later
saveRDS(fast_markers, "markers_raw.rds")
cat(sprintf("[%s] INFO | Raw markers saved: %d rows across %d clusters\n",
            Sys.time(), nrow(fast_markers), n_distinct(fast_markers$cluster)),
    file = stderr())

# ── Stage 2: Post-filter for enrichR input ───────────────────────────────────
# p_val_adj is essential here — return.thresh above used nominal p only
markers_clean <- fast_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 150, with_ties = FALSE) %>%
  ungroup()

# Save clean version — this is the golden standard RDS for enrichr_pipeline.R
saveRDS(markers_clean, "markers_clean.rds")
cat(sprintf("[%s] INFO | Clean markers saved: %d rows | median %.0f genes/cluster\n",
            Sys.time(), nrow(markers_clean),
            median(table(markers_clean$cluster))),
    file = stderr())

# ── Sanity check ─────────────────────────────────────────────────────────────
markers_clean %>%
  group_by(cluster) %>%
  summarise(
    n_genes    = n(),
    median_lfc = round(median(avg_log2FC), 2),
    min_padj   = signif(min(p_val_adj), 3)
  ) %>%
  print()

# Verify final_donor exists in Seurat metadata
cat("final_donor exists:", "Sample" %in% names(seurat_obj@meta.data), "\n")
cat("N unique donors:", length(unique(seurat_obj$Sample)), "\n")
head(seurat_obj$Sample)

counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

sce_clean <- SingleCellExperiment(
  assays  = list(counts = counts_mat),
  colData = DataFrame(
    Sample     = seurat_obj$Sample,
    cluster_id = seurat_obj$class,
    dummy      = "x"              # required — dreamlet 1.9 bug workaround
  )
)

pb <- dreamlet::aggregateToPseudoBulk(
  sce_clean,
  assay      = "counts",
  cluster_id = "cluster_id",
  sample_id  = "Sample",
  verbose    = TRUE
)

colData(pb)$Sample <- rownames(colData(pb))  # required — dreamlet 1.9 bug workaround

all_clusters <- assayNames(pb)

dreamlet_results <- map(all_clusters, function(ct) {
  cat(sprintf("Running: %s vs rest\n", ct))
  tryCatch({
    fit <- dreamletCompareClusters(pb, c(ct, "rest"), method = "fixed")
    df  <- topTable(fit, coef = "compare", number = Inf)
    df$gene    <- rownames(df)
    df$cluster <- ct
    rownames(df) <- NULL
    df
  }, error = function(e) {
    cat(sprintf("FAILED %s: %s\n", ct, e$message))
    NULL
  })
}) %>%
  set_names(all_clusters) %>%
  compact()

cat("Clusters completed:", paste(names(dreamlet_results), collapse=", "), "\n")
cat("Total DEGs:", nrow(bind_rows(dreamlet_results)), "\n")

# Save
dreamlet_merged <- bind_rows(dreamlet_results)
saveRDS(dreamlet_merged, "dreamlet_degs_all_clusters.rds")

# Quick sanity check on the DEGs before enrichR
dreamlet_merged %>%
  group_by(cluster) %>%
  summarise(
    total_genes   = n(),
    sig_up        = sum(adj.P.Val < 0.05 & logFC > 0.5),
    sig_down      = sum(adj.P.Val < 0.05 & logFC < -0.5)
  ) %>%
  print()

# Save and feed to enrichR pipeline
saveRDS(dreamlet_merged, "miBRAIN_dreamlet_degs_all_clusters.rds")

df_cts <- cellTypeSpecificity(pb)
# retain only genes with total CPM summed across cell type > 100
df_cts <- df_cts[df_cts$totalCPM > 100, ]
# Define N, the number of top genes you want (in this case, 5)
N <- 500

# Use apply() with MARGIN=2 (columns) and a custom function
# The function sorts the column and returns the top N row names (genes)
top_n_genes <- apply(df_cts, 2, function(col_values) {
    # 1. Use 'order' with 'decreasing=TRUE' to get the indices 
    #    of the values sorted from largest to smallest.
    # 2. Select the first N indices (the top N values).
    top_indices <- order(col_values, decreasing = TRUE)[1:N]
    
    # 3. Use these indices to select the corresponding row names (genes)
    #    from the full data frame/matrix's row names.
    return(rownames(df_cts)[top_indices])
})

# retain only genes with total CPM summed across cell type > 100
df_cts <- df_cts[df_cts$totalCPM > 100, ]

# Violin plot of specificity score for each cell type
plotViolin(df_cts)
genes <- rownames(df_cts)[apply(df_cts, 2, which.max)]
plotPercentBars(df_cts, genes = genes)
dreamlet::plotHeatmap(df_cts, genes = genes)

# Define N, the number of top genes you want (in this case, 5)
N <- 5

# Use apply() with MARGIN=2 (columns) and a custom function
# The function sorts the column and returns the top N row names (genes)
top_n_genes <- apply(df_cts, 2, function(col_values) {
    # 1. Use 'order' with 'decreasing=TRUE' to get the indices 
    #    of the values sorted from largest to smallest.
    # 2. Select the first N indices (the top N values).
    top_indices <- order(col_values, decreasing = TRUE)[1:N]
    
    # 3. Use these indices to select the corresponding row names (genes)
    #    from the full data frame/matrix's row names.
    return(rownames(df_cts)[top_indices])
})

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 200)
dreamlet::plotHeatmap(df_cts, genes = unique(as.vector(top_n_genes[,-1])))

head(seurat_obj@meta.data)

# Pull sample-level metadata from Seurat
sample_meta <- seurat_obj@meta.data %>%
  as.data.frame() %>%
  dplyr::select(Sample, dx,
                nCount_RNA, n_features) %>%
  distinct(Sample, .keep_all = TRUE)

rownames(sample_meta) <- sample_meta$Sample
sample_meta$Sample <- NULL

cat("sample_meta rows:", nrow(sample_meta), "\n")
print(head(sample_meta))

# Merge into pb colData
colData(pb) <- cbind(colData(pb), sample_meta[rownames(colData(pb)), ])

cat("colData columns after merge:\n")
print(colnames(colData(pb)))
print(head(colData(pb)))

form_varpart <- ~   (1|dx) +  n_features


# Normalize, filter, voom + weights
res_proc_class <- processAssays(
  pb,
  form_varpart,
  min.count = 5,BPPARAM = SnowParam(workers = 6, progressbar = TRUE)
)

vp_lst <- fitVarPart(res_proc_class, form_varpart,
                     BPPARAM = SnowParam(workers = 12, progressbar = TRUE))

plotVarPart(vp_lst, label.angle = 60)

# show voom plot for each cell clusters
plotVoom(res_proc_class)

# The specific indices or cluster names you want to retain
target_clusters <- c(1, 2,3,5,6,7)
# Use single brackets to extract a sub-list
res_proc_class_subset <- res_proc_class[target_clusters]

# show voom plot for each cell clusters
plotVoom(res_proc_class_subset)

CONTRASTS <- c(E3vsCntrl = "dxA53T_E3 - dxCntrl",E4vsCntrl = "dxA53T_E4 - dxCntrl",E4vsE3 = "dxA53T_E4 - dxA53T_E3",
              A53TvsCtrl = "(dxA53T_E4+dxA53T_E3)/2 - dxCntrl")

# --- Model 1: Baseline with main covariates only ---
form_dreamlet1 <- ~ n_features+dx + 0

res_dl_dx1 <- dreamlet(res_proc_class_subset, form_dreamlet1, 
                       contrasts = CONTRASTS,
                      BPPARAM = SnowParam(workers = 12, progressbar = TRUE))

library(dplyr)
library(qvalue)

summarize_de <- function(res, coef_name) {
  summary_tab <- topTable(res, coef = coef_name, number = Inf) %>%
    as_tibble() %>%
    group_by(assay) %>%
    summarize(
      nGenes = n(),
      
      # significance
      nDE = sum(adj.P.Val < 0.05, na.rm = TRUE),
      
      # directionality
      nUp = sum(adj.P.Val < 0.05 & logFC > 0, na.rm = TRUE),
      nDown = sum(adj.P.Val < 0.05 & logFC < 0, na.rm = TRUE),
      
      # effect-size aware counts (optional but better)
      nUp_strict = sum(adj.P.Val < 0.05 & logFC > 0.25, na.rm = TRUE),
      nDown_strict = sum(adj.P.Val < 0.05 & logFC < -0.25, na.rm = TRUE),
      
      # global signal
      pi1 = 1 - pi0est(P.Value)$pi0
    ) %>%
    mutate(assay = factor(assay, assayNames(res)))
  
  return(summary_tab)
}

res_E3  <- summarize_de(res_dl_dx1, "E3vsCntrl")
res_E4  <- summarize_de(res_dl_dx1, "E4vsCntrl")
res_E4vsE3 <- summarize_de(res_dl_dx1, "E4vsE3")

write_csv(res_E3,"Summary_DEGs_E3vsCntrl.csv")
write_csv(res_E4,"Summary_DEGs_E4vsCntrl.csv")
write_csv(res_E4vsE3,"Summary_DEGs_E4vsE3.csv")

res_E4

# --- Optional: Summarize DE Results ---
COEF = 'E3vsCntrl'
summary_tab <- topTable(res_dl_dx1, coef = COEF, number = Inf) %>%
  as_tibble() %>%
  group_by(assay) %>%
  summarize(
    nDE = sum(adj.P.Val < 0.05),
    pi1 = 1 - pi0est(P.Value)$pi0,
    nGenes = length(adj.P.Val)
  ) %>%
  mutate(assay = factor(assay, assayNames(res_dl_dx1)))
summary_tab

# --- Optional: Summarize DE Results ---
COEF = 'E4vsCntrl'
summary_tab <- topTable(res_dl_dx1, coef = COEF, number = Inf) %>%
  as_tibble() %>%
  group_by(assay) %>%
  summarize(
    nDE = sum(adj.P.Val < 0.05),
    pi1 = 1 - pi0est(P.Value)$pi0,
    nGenes = length(adj.P.Val)
  ) %>%
  mutate(assay = factor(assay, assayNames(res_dl_dx1)))
summary_tab

# --- Optional: Summarize DE Results ---
COEF = 'A53TvsCtrl'
summary_tab <- topTable(res_dl_dx1, coef = COEF, number = Inf) %>%
  as_tibble() %>%
  group_by(assay) %>%
  summarize(
    nDE = sum(adj.P.Val < 0.05),
    pi1 = 1 - pi0est(P.Value)$pi0,
    nGenes = length(adj.P.Val)
  ) %>%
  mutate(assay = factor(assay, assayNames(res_dl_dx1)))
summary_tab
# write.csv(summary_tab, "./summary_dxPD_CLASS_cova_DiffPDvsCTRL_raw.csv")

saveRDS(res_dl_dx1,'miBRAIN_CLASS_cova_Diff.rds')

# saveRDS(res_dl_dx1,'CLASS_cova_DiffPDvsCTRL_1.9.rds')

library(ggplot2)
library(dplyr)
library(ggrepel)
log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

#' Create a volcano plot for differential expression analysis.
#'
#' @param AllDF Data frame containing DE results.
#' @param cellType String specifying the cell type to filter.
#' @param n_genes Integer for the number of top UP/DOWN genes to label.
#' @param nominal Logical; if TRUE, uses 'P.Value' instead of 'adj.P.Val'.
#' @return A ggplot object representing the volcano plot.
my_Vlc_plot <- function(AllDF, cellType, n_genes = 10, nominal = FALSE,logfc_thres=0.5) {
  
  # 1. Pipeline Switch: Determine the target column and Y-axis label
  pval_col <- ifelse(nominal, "P.Value", "adj.P.Val")
  y_label  <- ifelse(nominal, "-log10(Nominal P-Value)", "-log10(Adjusted P-Value)")

  cat(sprintf("[%s] | INFO | Generating volcano plot for %s using %s\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cellType, pval_col), file=stderr())

  # 2. Filter data
  DF <- AllDF %>% dplyr::filter(assay == cellType)
  
  if (nrow(DF) == 0) {
    cat(sprintf("[%s] | WARN | No data found for cell type: %s. Returning NULL.\n", 
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cellType), file=stderr())
    return(NULL)
  }

  # Standardize the p-value column name internally so downstream code is clean
  DF <- DF %>% dplyr::rename(target_pval = all_of(pval_col))

  # 3. Categorize expression safely 
  DF <- DF %>%
    mutate(
      diffexpressed = case_when(
        logFC > logfc_thres & target_pval < 0.05 ~ "UP",
        logFC < -logfc_thres & target_pval < 0.05 ~ "DOWN",
        TRUE ~ "NS" 
      )
    )

  # 4. Safely label top genes (sort by the selected p-value)
  top_up <- DF %>% 
    filter(diffexpressed == "UP") %>% 
    arrange(target_pval) %>% 
    head(n_genes) %>% 
    pull(gene.name)
    
  top_down <- DF %>% 
    filter(diffexpressed == "DOWN") %>% 
    arrange(target_pval) %>% 
    head(n_genes) %>% 
    pull(gene.name)

  DF <- DF %>%
    mutate(delabel = ifelse(gene.name %in% c(top_up, top_down), gene.name, NA))

  # 5. Dynamic Limits
  x_limit <- ceiling(max(abs(DF$logFC), na.rm = TRUE))

  # 6. Plotting
  plot <- ggplot(DF, aes(x = logFC, y = -log10(target_pval), col = diffexpressed, label = delabel)) +
    geom_point(size = 0.8, alpha = 0.7) + 
    geom_text_repel(color = 'black', size = 3.5, max.overlaps = 15) + 
    scale_color_manual(values = c("UP" = "#C0392B", "DOWN" = "#2980B9", "NS" = "grey80")) + 
    geom_vline(xintercept = c(-logfc_thres, logfc_thres), col = "black", linetype = "dashed", alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed", alpha = 0.5) + 
    xlim(-x_limit, x_limit) + 
    theme_light() + 
    labs(
      title = cellType, 
      subtitle = paste0("Total genes tested: ", nrow(DF)),
      x = "log2(Fold Change)",
      y = y_label
    ) + 
    theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom", 
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)))

  cat(sprintf("[%s] | INFO | Plot generation complete for %s\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cellType), file=stderr())
              
  return(plot) 
}


library(dreamlet)
library(ggplot2)
library(tidyverse)
library(aplot)
library(ggtree)
library(scattermore)
library(zenith)
library(crumblr)
library(GSEABase)
library(qvalue)
library(BiocParallel)
library(cowplot)
library(presto)
library(DelayedArray)

# --- Print Versions for Reproducibility ---

cat("Session package versions:\n")
cat(paste0("dreamlet v", packageVersion("dreamlet")), "\n")
cat(paste0("crumblr v", packageVersion("crumblr")), "\n")
cat(paste0("variancePartition v", packageVersion("variancePartition")), "\n")
cat(paste0("zenith v", packageVersion("zenith")), "\n")
cat(paste0("zellkonverter v", packageVersion("zellkonverter")), "\n")
cat(paste0("BiocManager v", BiocManager::version()), "\n")

print(getwd())

# --- Step 1: Load H5AD and Metadata ---
res.dl <- readRDS('miBRAIN_CLASS_cova_Diff.rds')

Alldf <- topTable(res.dl, coef = 'E3vsCntrl', number = Inf)

getwd()

# ---- Main function ------------------------------------------------------------
run_fgsea_analysis <- function(gmt_list,
                               df         = AllDF,
                               coef_name  = "E3vsCntrl", # NEW: Pass your COEF here
                               top_n      = TOP_N,
                               top_n_up   = TOP_N_UP,
                               top_n_down = TOP_N_DOWN,
                               nperm      = NPERM,
                               out_prefix = "fgsea") {

  all_gmt    <- do.call(c, gmt_list)
  single_gmt <- length(gmt_list) == 1
  
  # Ensure df is a standard tibble
  df_clean <- tibble::as_tibble(as.data.frame(df))
  clusters <- unique(df_clean$assay)
  
  if (single_gmt) {
       names(all_gmt) <- sub("^.*\\.[^_]*_", "", names(all_gmt))
  }
  
  # NEW: Extract UP and DOWN group names dynamically from the coefficient string
  groups <- unlist(strsplit(coef_name, "vs", fixed = TRUE))
  group_up   <- if (length(groups) >= 1) groups[1] else "UP"
  group_down <- if (length(groups) >= 2) groups[2] else "DOWN"
  
  # ---- Run fgsea per cluster --------------------------------------------------
  fgsea_results <- lapply(clusters, function(cl) {
    ranked <- df_clean %>%
      dplyr::filter(assay == cl) %>%
      dplyr::arrange(desc(t)) %>%
      { setNames(.$t, .$ID) }

    fgsea::fgsea(pathways = all_gmt, stats = ranked, minSize = 15, maxSize = 500, nPermSimple = nperm) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(cluster = cl)
  })

  fgsea_all_df <- dplyr::bind_rows(fgsea_results) %>%
    dplyr::select(-leadingEdge) %>%
    dplyr::arrange(cluster, padj)

  fgsea_up_df   <- fgsea_all_df %>% dplyr::filter(NES > 0)
  fgsea_down_df <- fgsea_all_df %>% dplyr::filter(NES < 0)

  # ---- Top N slices -----------------------------------------------------------
  top_up   <- fgsea_up_df   %>% dplyr::group_by(cluster) %>% dplyr::slice_max(abs(NES), n = top_n,      with_ties = FALSE) %>% dplyr::ungroup()
  top_down <- fgsea_down_df %>% dplyr::group_by(cluster) %>% dplyr::slice_max(abs(NES), n = top_n,      with_ties = FALSE) %>% dplyr::ungroup()
  top_up_c <- fgsea_up_df   %>% dplyr::group_by(cluster) %>% dplyr::slice_max(abs(NES), n = top_n_up,   with_ties = FALSE) %>% dplyr::ungroup()
  top_dn_c <- fgsea_down_df %>% dplyr::group_by(cluster) %>% dplyr::slice_max(abs(NES), n = top_n_down, with_ties = FALSE) %>% dplyr::ungroup()

  # ---- Tile builder -----------------------------------------------------------
  make_tile <- function(top_df, full_df, title, scale = "diverging") {

    row_order <- top_df %>%
      dplyr::arrange(cluster, desc(NES)) %>%
      dplyr::mutate(pathway = stringr::str_trunc(pathway, 50)) %>%
      dplyr::pull(pathway) %>% unique()

    plot_df <- tidyr::expand_grid(pathway = unique(top_df$pathway),
                                  cluster = unique(full_df$cluster)) %>%
      dplyr::left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj),
                       by = c("pathway", "cluster")) %>%
      dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                    sig     = !is.na(padj) & padj < 0.05,
                    pathway = factor(pathway, levels = rev(row_order)))

    lim <- max(abs(plot_df$NES), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) { 
      cat(sprintf("[%s] WARNING | Empty plot: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), title), file = stderr())
      return(NULL) 
    }

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cluster, y = pathway, fill = NES)) +
      ggplot2::geom_tile(ggplot2::aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
      ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 40, hjust = 1),
                     axis.text.y  = ggplot2::element_text(size = 7),
                     panel.grid   = ggplot2::element_blank(),
                     plot.title   = ggplot2::element_text(face = "bold"))

    if (scale == "diverging") {
      p <- p + ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                             midpoint = 0, limits = c(-lim, lim),
                                             na.value = "grey85", name = "NES")
    } else if (scale == "red") {
      p <- p + ggplot2::scale_fill_gradient(low = "#fddbc7", high = "#b2182b",
                                            na.value = "grey85", name = "NES")
    } else if (scale == "blue") {
      p <- p + ggplot2::scale_fill_gradient(low = "#2166ac", high = "#d1e5f0",
                                            na.value = "grey85", name = "NES")
    }
    return(p)
  }

  # ---- Build plots ------------------------------------------------------------
  # NEW: Titles now dynamically call group_up and group_down instead of hardcoded strings
  p_up   <- make_tile(top_up,   fgsea_up_df,   paste0(out_prefix, ": ", group_up, "-enriched"),      scale = "red")
  p_down <- make_tile(top_down, fgsea_down_df, paste0(out_prefix, ": ", group_down, "-enriched"),    scale = "blue")

  # Combined
  top_combined       <- dplyr::bind_rows(top_up_c, top_dn_c)
  row_order_combined <- top_combined %>%
    dplyr::arrange(desc(NES > 0), cluster, desc(NES)) %>%
    dplyr::mutate(pathway = stringr::str_trunc(pathway, 50)) %>%
    dplyr::pull(pathway) %>% unique()

  plot_df_combined <- tidyr::expand_grid(
      pathway = unique(top_combined$pathway),
      cluster = unique(fgsea_all_df$cluster)) %>%
    dplyr::left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
                     by = c("pathway", "cluster")) %>%
    dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                  sig     = !is.na(padj) & padj < 0.05,
                  pathway = factor(pathway, levels = rev(row_order_combined)))

  lim_c <- max(abs(plot_df_combined$NES), na.rm = TRUE)

  p_combined <- ggplot2::ggplot(plot_df_combined, ggplot2::aes(x = cluster, y = pathway, fill = NES)) +
    ggplot2::geom_tile(ggplot2::aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                  midpoint = 0, limits = c(-lim_c, lim_c),
                                  na.value = "grey85", name = "NES") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = paste0(out_prefix, "\n ", group_up, " vs ", group_down), x = NULL, y = NULL) + # NEW: Dynamic combined title
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid  = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(face = "bold"))

  # Overall by |NES|
  top_ov <- fgsea_all_df %>% dplyr::group_by(cluster) %>% dplyr::slice_max(abs(NES), n = top_n, with_ties = FALSE) %>% dplyr::ungroup()
  p_overall <- make_tile(top_ov, fgsea_all_df, paste0(out_prefix, " \n Overall top pathways (|NES|)"), scale = "diverging")

  # ---- Save -------------------------------------------------------------------
  readr::write_csv(fgsea_up_df,   paste0(out_prefix, "_UP_df.csv"))
  readr::write_csv(fgsea_down_df, paste0(out_prefix, "_DOWN_df.csv"))

  ggplot2::ggsave(paste0(out_prefix, "_tile_UP.pdf"),        p_up,        width = 9,  height = 7,  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(out_prefix, "_tile_DOWN.pdf"),      p_down,      width = 9,  height = 7,  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(out_prefix, "_tile_combined.pdf"), p_combined, width = 10, height = 10, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(out_prefix, "_tile_overall.pdf"),  p_overall,  width = 10, height = 8,  device = grDevices::cairo_pdf)

  # Log completion to stderr
  cat(sprintf("[%s] INFO | Done: %s (%s vs %s) | UP: %d | DOWN: %d rows\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), out_prefix, group_up, group_down, nrow(fgsea_up_df), nrow(fgsea_down_df)), 
      file = stderr())

  invisible(list(all = fgsea_all_df, up = fgsea_up_df, down = fgsea_down_df,
                 p_up = p_up, p_down = p_down, p_combined = p_combined, p_overall = p_overall))
}

getwd()

library(tidyverse)
library(fgsea)

# ---- Parameters ---------------------------------------------------------------
TOP_N_UP   <- 5
TOP_N_DOWN <- 0
TOP_N      <- 5
NPERM      <- 1000

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt  <- gmtPathways(KEGG_GMT)
go_bp_gmt <- gmtPathways(GOBP_GMT)
go_mf_gmt <- gmtPathways(GOMF_GMT)

unique(Alldf$assay)

COEF <- 'E3vsCntrl'
Alldf <- as.data.frame(topTable(res.dl, coef = COEF, number = Inf)) %>% 
filter(assay %in% c('Astrocytes','Early neurons','Neurons','Progenitors/OPC','Endothelial cells','Mural/stromal cells'))

# Assuming you want to pull the gene names from an existing AllDF object, 
# or you could use rownames(Alldf) if standard limma formatting applies.
Alldf$gene.name <- Alldf$ID 

# ---- Run separately per GMT ---------------------------------------------------
res_kegg  <- run_fgsea_analysis(gmt_list = list(KEGG  = kegg_gmt), df = Alldf, coef_name = COEF, out_prefix = "kegg")
res_go_bp <- run_fgsea_analysis(gmt_list = list(GO_BP = go_bp_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_bp")
res_go_mf <- run_fgsea_analysis(gmt_list = list(GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_mf")

# ---- Or run combined ----------------------------------------------------------
res_all   <- run_fgsea_analysis(gmt_list = list(KEGG = kegg_gmt, GO_BP = go_bp_gmt, GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "all_gmt")

options(repr.plot.width = 7, repr.plot.height = 7, repr.plot.res = 200)
res_kegg$p_combined
res_go_bp$p_combined
res_go_mf$p_combined

head(res_go_bp$all %>% filter(cluster=='Neurons',padj<0.01) %>% arrange(desc(-log10(padj))),40)

Neurons_GOBP_Sign = res_go_bp$all %>% filter(cluster=='Neurons',padj<0.05) %>% arrange(desc(NES))
write_csv(Neurons_GOBP_Sign,'E3vsCntrl_Neurons_GOBP_Enrichment.csv')
getwd()

COEF <- 'E4vsE3'
Alldf <- topTable(res.dl, coef = COEF, number = Inf)

# Assuming you want to pull the gene names from an existing AllDF object, 
# or you could use rownames(Alldf) if standard limma formatting applies.
Alldf$gene.name <- Alldf$ID 

# ---- Run separately per GMT ---------------------------------------------------
res_kegg  <- run_fgsea_analysis(gmt_list = list(KEGG  = kegg_gmt), df = Alldf, coef_name = COEF, out_prefix = "kegg")
res_go_bp <- run_fgsea_analysis(gmt_list = list(GO_BP = go_bp_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_bp")
res_go_mf <- run_fgsea_analysis(gmt_list = list(GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_mf")

# ---- Or run combined ----------------------------------------------------------
res_all   <- run_fgsea_analysis(gmt_list = list(KEGG = kegg_gmt, GO_BP = go_bp_gmt, GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "all_gmt")

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
res_kegg$p_combined
res_go_bp$p_combined
res_go_mf$p_combined

COEF <- 'A53TvsCtrl'
Alldf <- topTable(res.dl, coef = COEF, number = Inf)

# Assuming you want to pull the gene names from an existing AllDF object, 
# or you could use rownames(Alldf) if standard limma formatting applies.
Alldf$gene.name <- Alldf$ID 

# ---- Run separately per GMT ---------------------------------------------------
res_kegg  <- run_fgsea_analysis(gmt_list = list(KEGG  = kegg_gmt), df = Alldf, coef_name = COEF, out_prefix = "kegg")
res_go_bp <- run_fgsea_analysis(gmt_list = list(GO_BP = go_bp_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_bp")
res_go_mf <- run_fgsea_analysis(gmt_list = list(GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "go_mf")

# ---- Or run combined ----------------------------------------------------------
res_all   <- run_fgsea_analysis(gmt_list = list(KEGG = kegg_gmt, GO_BP = go_bp_gmt, GO_MF = go_mf_gmt), df = Alldf, coef_name = COEF, out_prefix = "all_gmt")

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 300)
res_go_bp$p_combined
res_kegg$p_combined
res_go_mf$p_combined

coefNames(res.dl)

library(tidyverse)
library(fgsea)

# ---- Parameters ---------------------------------------------------------------
TOP_N_UP   <- 5
TOP_N_DOWN <- 1
TOP_N      <- 5
NPERM      <- 1000

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt  <- gmtPathways(KEGG_GMT)
go_bp_gmt <- gmtPathways(GOBP_GMT)
go_mf_gmt <- gmtPathways(GOMF_GMT)

COEF <- 'E3vsCntrl'
ranked_ExN <- topTable(res.dl, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == "Excitatory neurons") %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

# Verify
head(names(ranked_ExN))  # should show ALG13, DOLPP1 etc now

res <- fgsea::fgsea(pathways = kegg_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

p <- res %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E3", "FALSE" = "Ctrl"),
                    name = NULL) +
  labs(title    = "KEGG: ExN \n E3 vs Ctrl",
       subtitle = "Top 20 by |NES| * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

file_name <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.png",COEF, 'KEGG')
file_name_pdf <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.pdf",COEF, 'KEGG')
ggsave(p,filename = file_name,width = 7,height =5 ,dpi=300)
ggsave(p,filename = file_name_pdf,width = 7,height =5)
options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 300)
p

res_gobp%>%filter(pathway%in%
                  c('GOBP_SYNAPSE_ORGANIZATION',
                   'GOBP_NEURON_APOPTOTIC_PROCESS',
                   'GOBP_NEURON_PROJECTION_GUIDANCE',
                   'GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY',
                   'GOBP_REGULATION_OF_ACTIN_FILAMENT_ORGANIZATION',
                   'GOBP_REGULATION_OF_ACTIN_FILAMENT_BASED_PROCESS',
                   'GOBP_NEURON_MIGRATION',
                   'GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION',
                   'GOBP_NEURON_DEATH',
                   'GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT')
                 )

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

# 1. Clean and Filter
p_data <- res_gobp %>%
  dplyr::filter(!is.na(padj)) %>%
  # Remove common prefixes (adjust regex based on your GMT names)
  dplyr::mutate(pathway = stringr::str_remove_all(pathway, "^(GOBP_|GO_BIOLOGICAL_PROCESS_)")) %>%
  dplyr::mutate(pathway = stringr::str_replace_all(pathway, "_", " ")) %>%
  # Selection Criteria: Filter by size and significance, then pick top NES
  dplyr::filter(size >= 100) %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::slice_max(abs(NES), n = 15, with_ties = FALSE) %>%
  dplyr::mutate(#pathway = stringr::str_trunc(pathway, 60),
                pathway = forcats::fct_reorder(pathway, NES))

# 2. Plot
p <- ggplot(p_data, aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # Significance markers
  geom_point(aes(x = NES + sign(NES) * 0.1), 
             shape = 8, size = 1.5, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "A53T", "FALSE" = "Ctrl"),
                    name = "Direction") +
  labs(title    = "GO BP: ExN | A53T vs Ctrl",
       subtitle = "Top 20 by |NES| (padj < 0.01)",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))+ theme(
    # Align title and subtitle to the right
    plot.title = element_text(hjust = 1, face = "bold"),
    plot.subtitle = element_text(hjust = 1)
  )

# 3. Save
file_base <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s", COEF, 'GOBP')
ggsave(p, filename = paste0(file_base, ".png"), width = 8, height = 6, dpi = 300)
ggsave(p, filename = paste0(file_base, ".pdf"), width = 8, height = 6)

p

COEF <- 'E3vsCntrl'
# GO BP
ranked_ExN <- topTable(res.dl, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == "Excitatory neurons") %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

res_gobp <- fgsea::fgsea(pathways = go_bp_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

p<-res_gobp %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "A53T", "FALSE" = "Ctrl"),
                    name = NULL) +
  labs(title    = "GO BP: ExN | A53T vs Ctrl",
       subtitle = "Top 20 by |NES|   * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

file_name <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.png",COEF, 'GOBP')
file_name_pdf <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.pdf",COEF, 'GOBP')
ggsave(p,filename = file_name,width = 7,height =5 ,dpi=300)
ggsave(p,filename = file_name_pdf,width = 7,height =5 )

p

# GO MF
res_gomf <- fgsea::fgsea(pathways = go_mf_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

p <- res_gomf %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E3", "FALSE" = "Ctrl"),
                    name = NULL) +
  labs(title    = "GO MF: ExN \n E3 vs Ctrl",
       subtitle = "Top 20 by |NES|   * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

file_name <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.png",COEF, 'GOMF')
file_name_pdf <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.pdf",COEF, 'GOMF')
ggsave(p,filename = file_name,width = 7,height =5 ,dpi=300)
ggsave(p,filename = file_name_pdf,width = 7,height =5)

coefNames(res.dl)

# # NEURONS E4 vs E3

COEF <- 'E4vsE3'
# GO BP
ranked_ExN <- topTable(res.dl, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == "Excitatory neurons") %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

res_gobp <- fgsea::fgsea(pathways = go_bp_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 300)

p<- res_gobp %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E4", "FALSE" = "E3"),
                    name = NULL) +
  labs(title    = "GO BP: ExN | E4 vs E3",
       subtitle = "Top 20 by |NES|   * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))


file_name <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.png",COEF, 'GOBP')
file_name_pdf <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.pdf",COEF, 'GOBP')
ggsave(p,filename = file_name,width = 7,height =5 ,dpi=300)
ggsave(p,filename = file_name_pdf,width = 7,height =5)
p

# GO MF
res_gomf <- fgsea::fgsea(pathways = go_mf_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)


p <- res_gomf %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E4", "FALSE" = "E3"),
                    name = NULL) +
  labs(title    = "GO MF: ExN \n E4 vs E3",
       subtitle = "Top 20 by |NES|   * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

file_name <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.png",COEF, 'GOMF')
file_name_pdf <- sprintf("figures/Enrichment/Enrichment_ExN_%s_%s.pdf",COEF, 'GOMF')
ggsave(p,filename = file_name,width = 7,height =5 ,dpi=300)
ggsave(p,filename = file_name_pdf,width = 7,height =5)
p

COEF <- 'E4vsE3'
cellType='Astrocytes'
prefix_name = 'Astro'
options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 300)


# GO BP
ranked_ExN <- topTable(res.dl, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == cellType) %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

res_gobp <- fgsea::fgsea(pathways = go_bp_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

res_gobp %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E4", "FALSE" = "E3"),
                    name = NULL) +
  labs(title    = "GO BP: Astro \n E4 vs E3",
       subtitle = "Top 20 by |NES| * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

# GO MF
res_gomf <- fgsea::fgsea(pathways = go_mf_gmt, stats = ranked_ExN,
                    minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
  tibble::as_tibble() %>%
  dplyr::select(-dplyr::any_of("leadingEdge")) %>%
  dplyr::arrange(padj)

res_gomf %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::slice_max(abs(NES), n = 20, with_ties = FALSE) %>%
  dplyr::mutate(pathway = stringr::str_trunc(pathway, 50),
                pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot(aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(data = ~ dplyr::filter(.x, padj < 0.05),
             aes(x = NES + sign(NES) * 0.08),
             shape = 8, size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    labels = c("TRUE" = "E4", "FALSE" = "E3"),
                    name = NULL) +
  labs(title    = "GO MF: Astro \n E4 vs E3",
       subtitle = "Top 20 by |NES|   * padj < 0.05",
       x = "NES", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

library(ggplot2)
library(dplyr)

plot_enrichment_classic <- function(fgsea_res,
                                    top_n    = 20,
                                    padj_cut = 0.05,
                                    title    = "GSEA Enrichment") {
  
  # Pull the combined result from the returned list
  df <- fgsea_res$all %>%
    dplyr::filter(!is.na(padj)) %>%
    # Take top N by |NES|, ensuring both directions are represented
    dplyr::slice_max(abs(NES), n = top_n, with_ties = FALSE) %>%
    dplyr::mutate(
      pathway   = stringr::str_trunc(pathway, 55),
      pathway   = stringr::str_wrap(pathway, width = 45),
      direction = ifelse(NES > 0, "Up", "Down"),
      sig_label = ifelse(padj < padj_cut, sprintf("p.adj=%.3f", padj), "n.s."),
      pathway   = forcats::fct_reorder(pathway, NES)
    )
  
  ggplot(df, aes(x = NES, y = pathway, fill = direction)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey30") +
    # Significance stars/dots on bars
    geom_point(
      data = df %>% dplyr::filter(padj < padj_cut),
      aes(x = NES + sign(NES) * 0.05),
      shape = 8, size = 2, color = "black", show.legend = FALSE
    ) +
    scale_fill_manual(values = c("Up" = "#b2182b", "Down" = "#2166ac"),
                      name = "Direction") +
    # Optionally color by -log10(padj) instead — swap fill to that
    labs(
      title    = title,
      subtitle = sprintf("Top %d pathways by |NES| — * padj < %.2f", top_n, padj_cut),
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y    = element_text(size = 8, lineheight = 0.85),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top",
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(color = "grey40", size = 9)
    )
}

# ---- Plot -------------------------------------------------------------------
p_ExN_kegg <- plot_enrichment_classic(
  fgsea_res = res_kegg_ExN,
  top_n     = 20,
  title     = paste0("KEGG — Excitatory Neurons | ", COEF)
)

print(p_ExN_kegg)

ggsave("kegg_ExcitatoryNeurons_barplot.pdf", p_ExN_kegg,
       width = 10, height = 8, device = cairo_pdf)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

df_cellCounts <- pb

# Extract metadata and cross-tabulate donors (rows) by cell class (columns)
counts_table <- table(seurat_obj$final_donor, seurat_obj$class)

# Convert directly to the required flat data frame format
df_cellCounts <- as.data.frame.matrix(counts_table)

# Run crumblr
cobj <- crumblr(df_cellCounts)

# Pull sample-level metadata from Seurat
sample_meta <- seurat_obj@meta.data %>%
  as.data.frame() %>%
  dplyr::select(final_donor, sex, sample_set, dx,
                log1p_total_counts, log1p_n_genes_by_counts,
                pct_counts_mito, pct_counts_ribo) %>%
  distinct(final_donor, .keep_all = TRUE)

rownames(sample_meta) <- sample_meta$final_donor
sample_meta$final_donor <- NULL

cat("sample_meta rows:", nrow(sample_meta), "\n")
print(head(sample_meta))

# Merge into pb colData
colData(pb) <- cbind(colData(pb), sample_meta[rownames(colData(pb)), ])

cat("colData columns after merge:\n")
print(colnames(colData(pb)))
print(head(colData(pb)))

info <- as.data.frame(colData(pb))
info

library(reformulas)
library(lme4)
library(variancePartition)

# 2. Update formula to exclude 'Sample' since N=1 per level. 
# Added sex and sample_set to capture more batch/biological variance.
form <- ~ (1 | dx) + (1 | sex) + (1 | sample_set) +log1p_total_counts

# 3. Fit model and plot
vp <- fitExtractVarPartModel(cobj, form, info)
fig.vp <- plotPercentBars(vp)

fig.vp

info$ind  <- info$Sample

library(ggplot2)

# Perform PCA
# use crumblr::standardize() to get values with
# approximately equal sampling variance,
# which is a key property for downstream PCA and clustering analysis.
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, info, by = "row.names")

# Plot PCA
#   color by Subject
#   shape by Stimulated vs unstimulated
ggplot(df_pca, aes(PC1, PC2, color = as.character(ind), shape = dx)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_color_discrete(name = "Subject") +
  xlab("PC1") +
  ylab("PC2")

heatmap(cobj$E)

# Explicitly set the reference level so the model calculates dxBD relative to Control
info$dx <- factor(info$dx, levels = c("Control", "BD"))

# Ensure batch variables are also factors
info$sex <- as.factor(info$sex)
info$sample_set <- as.factor(info$sample_set)

form <- ~ dx + sex + sample_set

fit <- dream(cobj, form, info)

fit <- eBayes(fit)

# Now 'dxBD' exists because Control is the established baseline
res <- topTable(fit, coef = "dxBD", number = Inf)
head(res)

hcl <- buildClusterTreeFromPB(pb)

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, hcl, coef = "dxBD")

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

options(repr.plot.width =14, repr.plot.height = 6, repr.plot.res = 300)
plotTreeTestBeta(res) +
  theme(legend.position = "bottom", legend.box = "vertical") |
  plotForest(res, hide = FALSE) |
  fig.vp

getwd()

AllDF = as.data.frame(topTable(readRDS('2026_03_06_Scanpy_020_2dpass_annotated_CLASS_cova_DiffPDvsCTRL_1.9_corrected.rds'), coef = 'Diff', number = Inf))
AllDF$gene.name= AllDF$ID

figure_path <- OUTDIR
plot_type = 'Volc'
cellType='Fibroblast'
nominal = TRUE
logfc_thres=0.5
plot.height  =5
plot.width =5
plot.res = 200
options(repr.plot.width = plot.height, repr.plot.height = plot.height, repr.plot.res = plot.res)
p <- my_Vlc_plot(AllDF, cellType=cellType, n_genes = 10, nominal = TRUE,logfc_thres=0.5)
fname <- file.path(figure_path,sprintf("%s_%s_corrected.png", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
fname <- file.path(figure_path,sprintf("%s_%s.pdf", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
log_info(sprintf("Saved: %s", fname))
p

figure_path <- OUTDIR
plot_type = 'Volc'
cellType='Astro'
nominal = TRUE
logfc_thres=0.5
plot.height  =5
plot.width =5
plot.res = 200
options(repr.plot.width = plot.height, repr.plot.height = plot.height, repr.plot.res = plot.res)
p <- my_Vlc_plot(AllDF, cellType=cellType, n_genes = 10, nominal = TRUE,logfc_thres=0.5)
fname <- file.path(figure_path,sprintf("%s_%s.png", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
fname <- file.path(figure_path,sprintf("%s_%s.pdf", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
log_info(sprintf("Saved: %s", fname))
p

figure_path <- OUTDIR
plot_type = 'Volc'
cellType='Neuron'
nominal = TRUE
logfc_thres=0.5
plot.height  =5
plot.width =5
plot.res = 200
options(repr.plot.width = plot.height, repr.plot.height = plot.height, repr.plot.res = plot.res)
p <- my_Vlc_plot(AllDF, cellType=cellType, n_genes = 10, nominal = TRUE,logfc_thres=0.5)
fname <- file.path(figure_path,sprintf("%s_%s.png", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
fname <- file.path(figure_path,sprintf("%s_%s.pdf", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
log_info(sprintf("Saved: %s", fname))
p

figure_path <- OUTDIR
plot_type = 'Volc'
cellType='IPC'
nominal = TRUE
logfc_thres=0.5
plot.height  =5
plot.width =5
plot.res = 200
options(repr.plot.width = plot.height, repr.plot.height = plot.height, repr.plot.res = plot.res)
p <- my_Vlc_plot(AllDF, cellType=cellType, n_genes = 10, nominal = TRUE,logfc_thres=0.5)
fname <- file.path(figure_path,sprintf("%s_%s.png", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
fname <- file.path(figure_path,sprintf("%s_%s.pdf", plot_type, cellType, nominal))
ggsave(fname, plot = p, width = plot.width, height = plot.height,
       dpi = plot.res, units = "in")
log_info(sprintf("Saved: %s", fname))
p

library(tidyverse)
library(fgsea)

# ---- Parameters ---------------------------------------------------------------
TOP_N_UP   <- 5
TOP_N_DOWN <- 1
TOP_N      <- 5
NPERM      <- 1000

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt  <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_bp_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
go_mf_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.mf.v7.5.1.symbols.gmt")

# ---- Main function ------------------------------------------------------------
# gmt_list: a named list e.g. list(KEGG = kegg_gmt, GO_BP = go_bp_gmt)
#           or a single gmt   e.g. list(GO_MF = go_mf_gmt)
# If single GMT passed, pathway prefixes (KEGG_, GOBP_, GOMF_ etc) are stripped.

run_fgsea_analysis <- function(gmt_list,
                               df         = AllDF,
                               top_n      = TOP_N,
                               top_n_up   = TOP_N_UP,
                               top_n_down = TOP_N_DOWN,
                               nperm      = NPERM,
                               out_prefix = "fgsea") {

  all_gmt    <- do.call(c, gmt_list)
  single_gmt <- length(gmt_list) == 1
  clusters   <- unique(df$assay)

  # ---- Run fgsea per cluster --------------------------------------------------
  fgsea_results <- lapply(clusters, function(cl) {
    ranked <- df %>%
      filter(assay == cl) %>%
      arrange(desc(t)) %>%
      { setNames(.$t, .$ID) }

    fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = nperm) %>%
      as_tibble() %>%
      mutate(cluster = cl)
  })

  fgsea_all_df <- bind_rows(fgsea_results) %>%
    dplyr::select(-leadingEdge) %>%
    arrange(cluster, padj)

  if (single_gmt) {
    fgsea_all_df <- fgsea_all_df %>%
      mutate(pathway = pathway %>%
               str_remove("^KEGG_|^GOBP_|^GOMF_|^GOCC_|^REACTOME_|^WP_") %>%
               str_replace_all("_", " ") %>%
               str_to_sentence())
  }

  fgsea_up_df   <- fgsea_all_df %>% filter(NES > 0)
  fgsea_down_df <- fgsea_all_df %>% filter(NES < 0)

  # ---- Top N slices -----------------------------------------------------------
  top_up   <- fgsea_up_df   %>% group_by(cluster) %>% slice_max(abs(NES), n = top_n,      with_ties = FALSE) %>% ungroup()
  top_down <- fgsea_down_df %>% group_by(cluster) %>% slice_max(abs(NES), n = top_n,      with_ties = FALSE) %>% ungroup()
  top_up_c <- fgsea_up_df   %>% group_by(cluster) %>% slice_max(abs(NES), n = top_n_up,   with_ties = FALSE) %>% ungroup()
  top_dn_c <- fgsea_down_df %>% group_by(cluster) %>% slice_max(abs(NES), n = top_n_down, with_ties = FALSE) %>% ungroup()

  # ---- Tile builder -----------------------------------------------------------
  make_tile <- function(top_df, full_df, title, scale = "diverging") {

    row_order <- top_df %>%
      arrange(cluster, desc(NES)) %>%
      mutate(pathway = str_trunc(pathway, 50)) %>%
      pull(pathway) %>% unique()

    plot_df <- expand_grid(pathway = unique(top_df$pathway),
                           cluster = unique(full_df$cluster)) %>%
      left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj),
                by = c("pathway", "cluster")) %>%
      mutate(pathway = str_trunc(pathway, 50),
             sig     = !is.na(padj) & padj < 0.05,
             pathway = factor(pathway, levels = rev(row_order)))

    lim <- max(abs(plot_df$NES), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) { message("Empty plot: ", title); return(NULL) }

    p <- ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
      geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
      coord_fixed() +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x  = element_text(angle = 40, hjust = 1),
            axis.text.y  = element_text(size = 7),
            panel.grid   = element_blank(),
            plot.title   = element_text(face = "bold"))

    if (scale == "diverging") {
      p <- p + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                                    midpoint = 0, limits = c(-lim, lim),
                                    na.value = "grey85", name = "NES")
    } else if (scale == "red") {
      p <- p + scale_fill_gradient(low = "#fddbc7", high = "#b2182b",
                                   na.value = "grey85", name = "NES")
    } else if (scale == "blue") {
      p <- p + scale_fill_gradient(low = "#2166ac", high = "#d1e5f0",
                                   na.value = "grey85", name = "NES")
    }
    return(p)
  }

  # ---- Build plots ------------------------------------------------------------
  p_up   <- make_tile(top_up,   fgsea_up_df,   paste0(out_prefix, " | BD-enriched (NES > 0)"),      scale = "red")
  p_down <- make_tile(top_down, fgsea_down_df, paste0(out_prefix, " | Control-enriched (NES < 0)"), scale = "blue")

  # Combined
  top_combined       <- bind_rows(top_up_c, top_dn_c)
  row_order_combined <- top_combined %>%
    arrange(desc(NES > 0), cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>% unique()

  plot_df_combined <- expand_grid(
      pathway = unique(top_combined$pathway),
      cluster = unique(fgsea_all_df$cluster)) %>%
    left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
              by = c("pathway", "cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order_combined)))

  lim_c <- max(abs(plot_df_combined$NES), na.rm = TRUE)

  p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim_c, lim_c),
                         na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = paste0(out_prefix, " | BD vs Control — enriched pathways"), x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))

  # Overall by |NES|
  top_ov <- fgsea_all_df %>% group_by(cluster) %>% slice_max(abs(NES), n = top_n, with_ties = FALSE) %>% ungroup()
  p_overall <- make_tile(top_ov, fgsea_all_df, paste0(out_prefix, " | Overall top pathways (|NES|)"), scale = "diverging")

  # ---- Save -------------------------------------------------------------------
  write_csv(fgsea_up_df,   paste0(out_prefix, "_UP_df.csv"))
  write_csv(fgsea_down_df, paste0(out_prefix, "_DOWN_df.csv"))

  ggsave(paste0(out_prefix, "_tile_UP.pdf"),       p_up,       width = 9,  height = 7,  device = cairo_pdf)
  ggsave(paste0(out_prefix, "_tile_DOWN.pdf"),     p_down,     width = 9,  height = 7,  device = cairo_pdf)
  ggsave(paste0(out_prefix, "_tile_combined.pdf"), p_combined, width = 10, height = 10, device = cairo_pdf)
  ggsave(paste0(out_prefix, "_tile_overall.pdf"),  p_overall,  width = 10, height = 8,  device = cairo_pdf)

  cat(sprintf("Done: %s | UP: %d | DOWN: %d rows\n", out_prefix, nrow(fgsea_up_df), nrow(fgsea_down_df)))

  invisible(list(all = fgsea_all_df, up = fgsea_up_df, down = fgsea_down_df,
                 p_up = p_up, p_down = p_down, p_combined = p_combined, p_overall = p_overall))
}

library(tidyverse)
library(fgsea)

# ---- Parameters ---------------------------------------------------------------
TOP_N_UP   <- 5
TOP_N_DOWN <- 1
TOP_N      <- 5
NPERM      <- 1000

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt  <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_bp_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
go_mf_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.mf.v7.5.1.symbols.gmt")

# ---- Run separately per GMT ---------------------------------------------------
res_kegg  <- run_fgsea_analysis(list(KEGG  = kegg_gmt),  out_prefix = "kegg")
res_go_bp <- run_fgsea_analysis(list(GO_BP = go_bp_gmt), out_prefix = "go_bp")
res_go_mf <- run_fgsea_analysis(list(GO_MF = go_mf_gmt), out_prefix = "go_mf")

# ---- Or run combined ----------------------------------------------------------
res_all   <- run_fgsea_analysis(list(KEGG = kegg_gmt, GO_BP = go_bp_gmt, GO_MF = go_mf_gmt), out_prefix = "all_gmt")

# ---- View plots ---------------------------------------------------------------
res_go_bp$p_up
go_mf_gmt$p_down
res_kegg$p_combined

# Load Gene Ontology database
# use gene 'SYMBOL', or 'ENSEMBL' id
# use get_MSigDB() to load MSigDB
# Use Cellular Component (i.e. CC) to reduce run time here
go.gs <- get_GeneOntology("BP", to = "SYMBOL")

# Run zenith gene set analysis on result of dreamlet
res_zenith_bp <- zenith_gsa(res.dl, coef = "Diff", go.gs)

# examine results for each ell type and gene set
head(res_zenith_bp)

plotZenithResults(res_zenith_bp, 5, 0)

res.dl

# Load Gene Ontology database
# use get_MSigDB() to load MSigDB
go.gs_mf <- get_GeneOntology("MF", to = "SYMBOL")
# Run zenith gene set analysis on result of dreamlet
res_zenith_mf <- zenith_gsa(res.dl, coef = "Diff", go.gs_mf)
# examine results for each ell type and gene set
plotZenithResults(res_zenith_mf, 5, 0)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
my_Vlc_plot(AllDF, cellType='Astro', n_genes = 10, nominal = TRUE,logfc_thres=0.5)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
my_Vlc_plot(AllDF, cellType='Neuron', n_genes = 10, nominal = TRUE,logfc_thres=0.5)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
my_Vlc_plot(AllDF, cellType='IPC', n_genes = 10, nominal = TRUE,logfc_thres=0.5)

head(AllDF)

library(fgsea)
nominal=TRUE
sig_theshold =0.05

# ---- 6. Tile plot — color = NES (diverging), coord_fixed for equal tiles ------
make_fgsea_tile_plot <- function(enrich_df,
                                 top_n    = TOP_N,
                                 title    = "FGSEA Enrichment Tile Plot",
                                 padj_max = PADJ_CUTOFF) {

  if (nrow(enrich_df) == 0) { message("Empty DF — skipping."); return(NULL) }

  # Top N per cluster by |NES| among significant terms
  top_terms <- enrich_df %>%
    filter(padj <= padj_max) %>%
    group_by(cluster) %>%
    slice_max(abs(NES), n = top_n, with_ties = FALSE) %>%
    ungroup()

  if (nrow(top_terms) == 0) {
    message(sprintf("Relaxing threshold — showing top %d by |NES|.", top_n))
    top_terms <- enrich_df %>%
      group_by(cluster) %>%
      slice_max(abs(NES), n = top_n, with_ties = FALSE) %>%
      ungroup()
  }

  all_terms    <- unique(top_terms$Description)
  all_clusters <- unique(enrich_df$cluster)

  plot_df <- expand_grid(Description = all_terms, cluster = all_clusters) %>%
    left_join(
      enrich_df %>% select(Description, cluster, NES, padj, size, database),
      by = c("Description", "cluster")
    ) %>%
    left_join(
      top_terms %>% distinct(Description, database),
      by = "Description", suffix = c("", ".src")
    ) %>%
    mutate(
      Description  = str_trunc(Description, 52, ellipsis = "..."),
      sig          = !is.na(padj) & padj <= padj_max,
      database.src = replace_na(database.src, "Other")
    )

  nes_lim <- max(abs(plot_df$NES), na.rm = TRUE)
  nes_lim <- if (is.finite(nes_lim) && nes_lim > 0) nes_lim else 2

  ggplot(plot_df, aes(x = cluster, y = Description, fill = NES)) +

    geom_tile(aes(color = sig), linewidth = 0.45, width = 0.92, height = 0.92) +
    scale_color_manual(values = c("TRUE" = "grey15", "FALSE" = "white"), guide = "none") +

    scale_fill_gradient2(
      low      = "#2166ac",   # control-enriched
      mid      = "white",
      high     = "#b2182b",   # BD-enriched
      midpoint = 0,
      limits   = c(-nes_lim, nes_lim),
      na.value = "grey88",
      name     = "NES"
    ) +

    # Separate facet strip per database — preserves visual grouping
    facet_grid(rows = vars(database.src), scales = "free_y", space = "free_y") +

    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +

    # Square tiles: 1 data unit wide = 1 data unit tall
    coord_fixed(ratio = 1) +

    labs(
      title    = title,
      subtitle = sprintf(
        "Top %d terms per cluster | red = BD-enriched (NES > 0) | blue = control-enriched (NES < 0)\nDark border = padj < %.2f | grey = absent in cluster",
        top_n, padj_max),
      x = "Cell cluster / assayID",
      y = NULL
    ) +

    theme_minimal(base_size = 11) +
    theme(
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(size = 8.5, color = "grey35", lineheight = 1.3),
      strip.text.y      = element_text(angle = 0, hjust = 0, size = 8,
                                       face = "bold", color = "grey20"),
      strip.background  = element_rect(fill = "grey92", color = NA),
      axis.text.x       = element_text(angle = 40, hjust = 1, size = 9),
      axis.text.y       = element_text(size = 7.5),
      panel.grid        = element_blank(),
      panel.spacing.y   = unit(0.35, "lines"),
      legend.position   = "right",
      legend.key.height = unit(1.4, "cm"),
      legend.title      = element_text(size = 9)
    )
}

universe <- unique(AllDF$ID)

if (nominal){    
    Up_Regulated <- AllDF %>% filter(logFC>0, P.Value<sig_theshold)
    Down_Regulated <- AllDF %>% filter(logFC<0, P.Value<sig_theshold)
}else{
    Up_Regulated <- AllDF %>% filter(logFC>0, adj.P.Val<sig_theshold)
    Down_Regulated <- AllDF %>% filter(logFC<0, adj.P.Val<sig_theshold)
}

bd_up_split_list <- split(Up_Regulated$ID,Up_Regulated$assay)
bd_down_split_list <- split(Down_Regulated$ID,Down_Regulated$assay)

# ============================================================
#  FGSEA Enrichment — flat script, no functions
#  Continuing from user's filtering + split code
# ============================================================

library(tidyverse)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(patchwork)

# ---- Parameters ----------------------------------------------------------------
sig_threshold <- 0.05    # p-value / adj.P.Val cutoff for ORA gene lists
nominal       <- TRUE   # TRUE = use P.Value; FALSE = use adj.P.Val
TOP_N         <- 5       # top N terms per cluster in tile plot
MIN_SIZE      <- 5
MAX_SIZE      <- 500
NPERM         <- 1000
SPECIES       <- "Homo sapiens"

# ---- Your existing code (bugs fixed) ------------------------------------------
universe <- unique(AllDF$ID)

if (nominal) {
  Up_Regulated   <- AllDF %>% filter(logFC > 0, P.Value   < sig_threshold)
  Down_Regulated <- AllDF %>% filter(logFC < 0, P.Value   < sig_threshold)  # BUG FIX: was logFC > 0
} else {
  Up_Regulated   <- AllDF %>% filter(logFC > 0, adj.P.Val < sig_threshold)
  Down_Regulated <- AllDF %>% filter(logFC < 0, adj.P.Val < sig_threshold)  # BUG FIX: was logFC > 0
}

bd_up_split_list   <- split(Up_Regulated$ID,   Up_Regulated$assay)
bd_down_split_list <- split(Down_Regulated$ID, Down_Regulated$assay)

# ---- Build ranked lists per cluster (full AllDF, ranked by Bz.std) ------------
# fgsea needs ALL genes ranked — not just the significant ones.
# bd_up/down_split_list defined above will be used to cross-check leading edges.

clusters <- unique(AllDF$assay)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- parse_gmt("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.entrez.gmt")
go_gmt   <- parse_gmt("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.entrez.gmt")
all_gmt  <- c(kegg_gmt, go_gmt)

colnames(AllDF)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})
fgsea_all_df   <- bind_rows(fgsea_results)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

library(tidyverse)
library(fgsea)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
all_gmt  <- c(kegg_gmt, go_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {
  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05)

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched \n(NES > 0)")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched \n(NES < 0)")

# ---- Save --------------------------------------------------------------------
# write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
# write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

# ggsave("tile_UP.pdf",   p_up,   width = 9, height = 7, device = cairo_pdf)
# ggsave("tile_DOWN.pdf", p_down, width = 9, height = 7, device = cairo_pdf)

# ---- Debug: check fgsea output ------------------------------------------------
cat("fgsea_all_df rows:", nrow(fgsea_all_df), "\n")
cat("Clusters found:", paste(unique(fgsea_all_df$cluster), collapse=", "), "\n")
cat("NES range:", range(fgsea_all_df$NES, na.rm=TRUE), "\n")
cat("UP rows:", nrow(fgsea_up_df), "\n")
cat("DOWN rows:", nrow(fgsea_down_df), "\n")
cat("Sample ranked vector (first cluster):\n")
cl1 <- unique(AllDF$assay)[1]
ranked_check <- AllDF %>% filter(assay == cl1) %>% arrange(desc(t)) %>% { setNames(.$t, .$ID) }
cat("  Length:", length(ranked_check), "\n")
cat("  First 3:", paste(head(names(ranked_check), 3), collapse=", "), "\n")
cat("  GMT first pathway genes (first 3):", paste(head(all_gmt[[1]], 3), collapse=", "), "\n")
cat("  Overlap with first pathway:", length(intersect(names(ranked_check), all_gmt[[1]])), "\n")

p_down

p_up

library(tidyverse)
library(fgsea)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
all_gmt  <- c(kegg_gmt, go_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched \n(NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Cntrl-enriched\n(NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
top_combined <- bind_rows(top_up, top_down)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control — all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 200)
p_up

all_gmt    <- do.call(c, gmt_list)
single_gmt <- length(gmt_list) == 1
clusters   <- unique(df$assay)

# ---- Run fgsea per cluster --------------------------------------------------
fgsea_results <- lapply(clusters, function(cl) {
ranked <- df %>%
  filter(assay == cl) %>%
  arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = nperm) %>%
  as_tibble() %>%
  mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>%
dplyr::select(-leadingEdge) %>%
arrange(cluster, padj)

if (single_gmt) {
fgsea_all_df <- fgsea_all_df %>%
  mutate(pathway = pathway %>%
           str_remove("^KEGG_|^GOBP_|^GOMF_|^GOCC_|^REACTOME_|^WP_") %>%
           str_replace_all("_", " ") %>%
           str_to_sentence())
}

fgsea_up_df   <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df <- fgsea_all_df %>% filter(NES < 0)

library(tidyverse)
library(fgsea)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
all_gmt  <- c(kegg_gmt, go_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched (NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched (NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
top_combined <- bind_rows(top_up, top_down)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control — all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways per cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)


p_up
p_down
p_combined
p_overall

library(tidyverse)
library(fgsea)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
all_gmt  <- c(kegg_gmt, go_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched (NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched (NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
top_combined <- bind_rows(top_up, top_down)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control — all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways per cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)

# ---- ORA on all DE genes (UP + DOWN merged) per cluster ----------------------
# Question: what pathways are dysregulated overall in BD, ignoring direction?
# Uses fgsea in "ORA mode" via fora() which takes a gene set + universe

all_DE <- bind_rows(Up_Regulated, Down_Regulated)

clusters <- unique(all_DE$assay)

fora_results <- lapply(clusters, function(cl) {
  genes_cl    <- all_DE %>% filter(assay == cl) %>% pull(ID) %>% unique()
  universe_cl <- AllDF %>% filter(assay == cl) %>% pull(ID) %>% unique()

  fora(pathways = all_gmt,
       genes    = genes_cl,
       universe = universe_cl,
       minSize  = 5) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fora_all_df <- bind_rows(fora_results) %>%
  arrange(cluster, padj)

# Top N per cluster by overlap (pval)
top_fora <- fora_all_df %>%
  group_by(cluster) %>%
  slice_min(padj, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_fora <- top_fora %>%
  arrange(cluster, padj) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_fora <- expand_grid(
    pathway = unique(top_fora$pathway),
    cluster = unique(fora_all_df$cluster)
  ) %>%
  left_join(fora_all_df %>% dplyr::select(pathway, cluster, padj, overlap),
            by = c("pathway", "cluster")) %>%
  mutate(pathway   = str_trunc(pathway, 50),
         sig       = !is.na(padj) & padj < 0.05,
         neg_log10 = -log10(padj + 1e-300),
         pathway   = factor(pathway, levels = rev(row_order_fora)))

p_fora <- ggplot(plot_df_fora, aes(x = cluster, y = pathway, fill = neg_log10)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient(low = "white", high = "#b2182b",
                      na.value = "grey85", name = "-log10(padj)") +
  coord_fixed() +
  labs(title = "Overall dysregulated pathways (UP + DOWN DE genes combined)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

write_csv(fora_all_df, "fora_overall_df.csv")
ggsave("tile_fora_overall.pdf", p_fora, width = 10, height = 8, device = cairo_pdf)


p_up
p_down
p_combined
p_overall
p_fora

library(tidyverse)
library(fgsea)

# ---- Load GMT files -----------------------------------------------------------
kegg_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")
# go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
all_gmt  <- c(kegg_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N <- 5

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched (NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched (NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
top_combined <- bind_rows(top_up, top_down)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control — all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways per cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)

# ---- ORA on all DE genes (UP + DOWN merged) per cluster ----------------------
# Question: what pathways are dysregulated overall in BD, ignoring direction?
# Uses fgsea in "ORA mode" via fora() which takes a gene set + universe

all_DE <- bind_rows(Up_Regulated, Down_Regulated)

clusters <- unique(all_DE$assay)

fora_results <- lapply(clusters, function(cl) {
  genes_cl    <- all_DE %>% filter(assay == cl) %>% pull(ID) %>% unique()
  universe_cl <- AllDF %>% filter(assay == cl) %>% pull(ID) %>% unique()

  fora(pathways = all_gmt,
       genes    = genes_cl,
       universe = universe_cl,
       minSize  = 5) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fora_all_df <- bind_rows(fora_results) %>%
  arrange(cluster, padj)

# Top N per cluster by overlap (pval)
top_fora <- fora_all_df %>%
  group_by(cluster) %>%
  slice_min(padj, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_fora <- top_fora %>%
  arrange(cluster, padj) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_fora <- expand_grid(
    pathway = unique(top_fora$pathway),
    cluster = unique(fora_all_df$cluster)
  ) %>%
  left_join(fora_all_df %>% dplyr::select(pathway, cluster, padj, overlap),
            by = c("pathway", "cluster")) %>%
  mutate(pathway   = str_trunc(pathway, 50),
         sig       = !is.na(padj) & padj < 0.05,
         neg_log10 = -log10(padj + 1e-300),
         pathway   = factor(pathway, levels = rev(row_order_fora)))

p_fora <- ggplot(plot_df_fora, aes(x = cluster, y = pathway, fill = neg_log10)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient(low = "white", high = "#b2182b",
                      na.value = "grey85", name = "-log10(padj)") +
  coord_fixed() +
  labs(title = "Overall dysregulated pathways (UP + DOWN DE genes combined)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

write_csv(fora_all_df, "fora_overall_df.csv")
ggsave("tile_fora_overall.pdf", p_fora, width = 10, height = 8, device = cairo_pdf)


p_up
p_down
p_combined
p_overall
p_fora

# library(tidyverse)
# library(fgsea)

# # ---- Load GMT files -----------------------------------------------------------
# MF_gmt <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.mf.v7.5.1.symbols.gmt")
# # go_gmt   <- gmtPathways("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
# all_gmt  <- c(kegg_gmt)

# # ---- Run fgsea per cluster ----------------------------------------------------
# clusters <- unique(AllDF$assay)

# fgsea_results <- lapply(clusters, function(cl) {
#   ranked <- AllDF %>%
#     filter(assay == cl) %>%
#     arrange(desc(t)) %>%
#     { setNames(.$t, .$ID) }

#   fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
#     as_tibble() %>%
#     mutate(cluster = cl)
# })

# fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
# fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
# fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# # ---- Top N per cluster --------------------------------------------------------
# TOP_N <- 5

# top_up <- fgsea_up_df %>%
#   group_by(cluster) %>%
#   slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
#   ungroup()

# top_down <- fgsea_down_df %>%
#   group_by(cluster) %>%
#   slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
#   ungroup()

# # ---- Tile plot ----------------------------------------------------------------
# make_tile <- function(top_df, full_df, title) {

#   # Order rows: group by nominating cluster, then sort by NES within group
#   row_order <- top_df %>%
#     arrange(cluster, desc(NES)) %>%
#     mutate(pathway = str_trunc(pathway, 50)) %>%
#     pull(pathway) %>%
#     unique()

#   plot_df <- expand_grid(pathway = unique(top_df$pathway),
#                          cluster = unique(full_df$cluster)) %>%
#     left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
#     mutate(pathway = str_trunc(pathway, 50),
#            sig     = !is.na(padj) & padj < 0.05,
#            pathway = factor(pathway, levels = rev(row_order)))

#   lim <- max(abs(plot_df$NES), na.rm = TRUE)

#   ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
#     geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
#     scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
#     scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
#                          midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
#     coord_fixed() +
#     labs(title = title, x = NULL, y = NULL) +
#     theme_minimal(base_size = 11) +
#     theme(axis.text.x = element_text(angle = 40, hjust = 1),
#           axis.text.y = element_text(size = 7),
#           panel.grid  = element_blank(),
#           plot.title  = element_text(face = "bold"))
# }

# p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched (NES > 0)") +
#   scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
# p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched (NES < 0)") +
#   scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# # ---- Combined UP + DOWN tile plot --------------------------------------------
# top_combined <- bind_rows(top_up, top_down)

# # Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
# row_order_combined <- top_combined %>%
#   arrange(desc(NES > 0), cluster, desc(NES)) %>%
#   mutate(pathway = str_trunc(pathway, 50)) %>%
#   pull(pathway) %>%
#   unique()

# plot_df_combined <- expand_grid(
#     pathway = unique(top_combined$pathway),
#     cluster = unique(fgsea_all_df$cluster)
#   ) %>%
#   left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
#             by = c("pathway", "cluster")) %>%
#   mutate(pathway = str_trunc(pathway, 50),
#          sig     = !is.na(padj) & padj < 0.05,
#          pathway = factor(pathway, levels = rev(row_order_combined)))

# lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

# p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
#   geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
#   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
#   scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
#                        midpoint = 0, limits = c(-lim_combined, lim_combined),
#                        na.value = "grey85", name = "NES") +
#   coord_fixed() +
#   labs(title = "BD vs Control — all enriched pathways", x = NULL, y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 40, hjust = 1),
#         axis.text.y = element_text(size = 7),
#         panel.grid  = element_blank(),
#         plot.title  = element_text(face = "bold"))

# # ---- Save --------------------------------------------------------------------
# write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
# write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

# ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
# ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
# ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# # ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# # fgsea_all_df already has results for all pathways in both directions.
# # Just pick top N by |NES| per cluster regardless of sign.

# top_overall <- fgsea_all_df %>%
#   group_by(cluster) %>%
#   slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
#   ungroup()

# row_order_overall <- top_overall %>%
#   arrange(cluster, desc(NES)) %>%
#   mutate(pathway = str_trunc(pathway, 50)) %>%
#   pull(pathway) %>%
#   unique()

# plot_df_overall <- expand_grid(
#     pathway = unique(top_overall$pathway),
#     cluster = unique(fgsea_all_df$cluster)
#   ) %>%
#   left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
#             by = c("pathway", "cluster")) %>%
#   mutate(pathway = str_trunc(pathway, 50),
#          sig     = !is.na(padj) & padj < 0.05,
#          pathway = factor(pathway, levels = rev(row_order_overall)))

# lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

# p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
#   geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
#   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
#   scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
#                        midpoint = 0, limits = c(-lim_overall, lim_overall),
#                        na.value = "grey85", name = "NES") +
#   coord_fixed() +
#   labs(title = "Overall most affected pathways per cluster (top 5 by |NES|)",
#        x = NULL, y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 40, hjust = 1),
#         axis.text.y = element_text(size = 7),
#         panel.grid  = element_blank(),
#         plot.title  = element_text(face = "bold"))

# ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)

# # ---- ORA on all DE genes (UP + DOWN merged) per cluster ----------------------
# # Question: what pathways are dysregulated overall in BD, ignoring direction?
# # Uses fgsea in "ORA mode" via fora() which takes a gene set + universe

# all_DE <- bind_rows(Up_Regulated, Down_Regulated)

# clusters <- unique(all_DE$assay)

# fora_results <- lapply(clusters, function(cl) {
#   genes_cl    <- all_DE %>% filter(assay == cl) %>% pull(ID) %>% unique()
#   universe_cl <- AllDF %>% filter(assay == cl) %>% pull(ID) %>% unique()

#   fora(pathways = all_gmt,
#        genes    = genes_cl,
#        universe = universe_cl,
#        minSize  = 5) %>%
#     as_tibble() %>%
#     mutate(cluster = cl)
# })

# fora_all_df <- bind_rows(fora_results) %>%
#   arrange(cluster, padj)

# # Top N per cluster by overlap (pval)
# top_fora <- fora_all_df %>%
#   group_by(cluster) %>%
#   slice_min(padj, n = TOP_N, with_ties = FALSE) %>%
#   ungroup()

# row_order_fora <- top_fora %>%
#   arrange(cluster, padj) %>%
#   mutate(pathway = str_trunc(pathway, 50)) %>%
#   pull(pathway) %>%
#   unique()

# plot_df_fora <- expand_grid(
#     pathway = unique(top_fora$pathway),
#     cluster = unique(fora_all_df$cluster)
#   ) %>%
#   left_join(fora_all_df %>% dplyr::select(pathway, cluster, padj, overlap),
#             by = c("pathway", "cluster")) %>%
#   mutate(pathway   = str_trunc(pathway, 50),
#          sig       = !is.na(padj) & padj < 0.05,
#          neg_log10 = -log10(padj + 1e-300),
#          pathway   = factor(pathway, levels = rev(row_order_fora)))

# p_fora <- ggplot(plot_df_fora, aes(x = cluster, y = pathway, fill = neg_log10)) +
#   geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
#   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
#   scale_fill_gradient(low = "white", high = "#b2182b",
#                       na.value = "grey85", name = "-log10(padj)") +
#   coord_fixed() +
#   labs(title = "Overall dysregulated pathways (UP + DOWN DE genes combined)",
#        x = NULL, y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 40, hjust = 1),
#         axis.text.y = element_text(size = 7),
#         panel.grid  = element_blank(),
#         plot.title  = element_text(face = "bold"))

# write_csv(fora_all_df, "fora_overall_df.csv")
# ggsave("tile_fora_overall.pdf", p_fora, width = 10, height = 8, device = cairo_pdf)


# p_up
# p_down
# p_combined
# p_overall
# p_fora

# # GO MF

library(tidyverse)
library(fgsea)


all_gmt  <- c(go_mf_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N_UP   <- 5   # terms per cluster in combined plot (BD-enriched)
TOP_N_DOWN <- 1   # terms per cluster in combined plot (control-enriched)
TOP_N      <- 5   # terms per cluster in UP/DOWN individual plots

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched\n(NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched\n(NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
# Use separate N for UP and DOWN in combined plot
top_up_combined   <- fgsea_up_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_UP,   with_ties = FALSE) %>% ungroup()
top_down_combined <- fgsea_down_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_DOWN, with_ties = FALSE) %>% ungroup()
top_combined <- bind_rows(top_up_combined, top_down_combined)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control \n all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways \nper cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)


# ---- FGSEA on merged UP+DOWN DE genes ranked by abs(t) -----------------------
# Question: what pathways are most dysregulated overall, ignoring direction?
# Rank by abs(t) — most dysregulated genes (either direction) at the top

all_DE <- bind_rows(Up_Regulated, Down_Regulated)

fgsea_overall_results <- lapply(clusters, function(cl) {
  ranked <- all_DE %>%
    filter(assay == cl) %>%
    arrange(desc(abs(t))) %>%
    { setNames(abs(.$t), .$ID) }

  fgsea(pathways    = all_gmt,
        stats       = ranked,
        minSize     = 5,
        maxSize     = 500,
        nPermSimple = NPERM) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_overall_df <- bind_rows(fgsea_overall_results) %>%
  dplyr::select(-leadingEdge) %>%
  arrange(cluster, padj)

# Top N per cluster by NES
top_overall2 <- fgsea_overall_df %>%
  group_by(cluster) %>%
  slice_max(NES, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall2 <- top_overall2 %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall2 <- expand_grid(
    pathway = unique(top_overall2$pathway),
    cluster = unique(fgsea_overall_df$cluster)
  ) %>%
  left_join(fgsea_overall_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall2)))

lim_overall2 <- max(abs(plot_df_overall2$NES), na.rm = TRUE)

p_overall2 <- ggplot(plot_df_overall2, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall2, lim_overall2),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall dysregulated pathways \n(ranked by |t|, UP+DOWN merged)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

write_csv(fgsea_overall_df, "fgsea_overall_df.csv")
ggsave("tile_overall.pdf", p_overall2, width = 10, height = 8, device = cairo_pdf)

p_up
p_down
p_combined
p_overall
p_overall2

# # GO BP

library(tidyverse)
library(fgsea)


all_gmt  <- c(go_bp_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N_UP   <- 5   # terms per cluster in combined plot (BD-enriched)
TOP_N_DOWN <- 1   # terms per cluster in combined plot (control-enriched)
TOP_N      <- 5   # terms per cluster in UP/DOWN individual plots

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched\n(NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched\n(NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
# Use separate N for UP and DOWN in combined plot
top_up_combined   <- fgsea_up_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_UP,   with_ties = FALSE) %>% ungroup()
top_down_combined <- fgsea_down_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_DOWN, with_ties = FALSE) %>% ungroup()
top_combined <- bind_rows(top_up_combined, top_down_combined)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control \n all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways \nper cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)


# ---- FGSEA on merged UP+DOWN DE genes ranked by abs(t) -----------------------
# Question: what pathways are most dysregulated overall, ignoring direction?
# Rank by abs(t) — most dysregulated genes (either direction) at the top

all_DE <- bind_rows(Up_Regulated, Down_Regulated)

fgsea_overall_results <- lapply(clusters, function(cl) {
  ranked <- all_DE %>%
    filter(assay == cl) %>%
    arrange(desc(abs(t))) %>%
    { setNames(abs(.$t), .$ID) }

  fgsea(pathways    = all_gmt,
        stats       = ranked,
        minSize     = 5,
        maxSize     = 500,
        nPermSimple = NPERM) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_overall_df <- bind_rows(fgsea_overall_results) %>%
  dplyr::select(-leadingEdge) %>%
  arrange(cluster, padj)

# Top N per cluster by NES
top_overall2 <- fgsea_overall_df %>%
  group_by(cluster) %>%
  slice_max(NES, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall2 <- top_overall2 %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall2 <- expand_grid(
    pathway = unique(top_overall2$pathway),
    cluster = unique(fgsea_overall_df$cluster)
  ) %>%
  left_join(fgsea_overall_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall2)))

lim_overall2 <- max(abs(plot_df_overall2$NES), na.rm = TRUE)

p_overall2 <- ggplot(plot_df_overall2, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall2, lim_overall2),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall dysregulated pathways \n(ranked by |t|, UP+DOWN merged)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

write_csv(fgsea_overall_df, "fgsea_overall_df.csv")
ggsave("tile_overall.pdf", p_overall2, width = 10, height = 8, device = cairo_pdf)

p_up
p_down
p_combined
p_overall
p_overall2

# # RUN KEGG

library(tidyverse)
library(fgsea)


all_gmt  <- c(kegg_gmt)

# ---- Run fgsea per cluster ----------------------------------------------------
clusters <- unique(AllDF$assay)

fgsea_results <- lapply(clusters, function(cl) {
  ranked <- AllDF %>%
    filter(assay == cl) %>%
    arrange(desc(t)) %>%
    { setNames(.$t, .$ID) }

  fgsea(pathways = all_gmt, stats = ranked, minSize = 5, maxSize = 500, nPermSimple = 1000) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_all_df <- bind_rows(fgsea_results) %>% dplyr::select(-leadingEdge)
fgsea_up_df    <- fgsea_all_df %>% filter(NES > 0)
fgsea_down_df  <- fgsea_all_df %>% filter(NES < 0)

# ---- Top N per cluster --------------------------------------------------------
TOP_N_UP   <- 5   # terms per cluster in combined plot (BD-enriched)
TOP_N_DOWN <- 1   # terms per cluster in combined plot (control-enriched)
TOP_N      <- 5   # terms per cluster in UP/DOWN individual plots

top_up <- fgsea_up_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

top_down <- fgsea_down_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- Tile plot ----------------------------------------------------------------
make_tile <- function(top_df, full_df, title) {

  # Order rows: group by nominating cluster, then sort by NES within group
  row_order <- top_df %>%
    arrange(cluster, desc(NES)) %>%
    mutate(pathway = str_trunc(pathway, 50)) %>%
    pull(pathway) %>%
    unique()

  plot_df <- expand_grid(pathway = unique(top_df$pathway),
                         cluster = unique(full_df$cluster)) %>%
    left_join(full_df %>% dplyr::select(pathway, cluster, NES, padj), by = c("pathway","cluster")) %>%
    mutate(pathway = str_trunc(pathway, 50),
           sig     = !is.na(padj) & padj < 0.05,
           pathway = factor(pathway, levels = rev(row_order)))

  lim <- max(abs(plot_df$NES), na.rm = TRUE)

  ggplot(plot_df, aes(x = cluster, y = pathway, fill = NES)) +
    geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = "NES") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold"))
}

p_up   <- make_tile(top_up,   fgsea_up_df,   "BD-enriched\n(NES > 0)") +
  scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey85", name = "NES")
p_down <- make_tile(top_down, fgsea_down_df, "Control-enriched\n(NES < 0)") +
  scale_fill_gradient(low = "#2166ac", high = "#d1e5f0", na.value = "grey85", name = "NES")



# ---- Combined UP + DOWN tile plot --------------------------------------------
# Use separate N for UP and DOWN in combined plot
top_up_combined   <- fgsea_up_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_UP,   with_ties = FALSE) %>% ungroup()
top_down_combined <- fgsea_down_df %>% group_by(cluster) %>% slice_max(abs(NES), n = TOP_N_DOWN, with_ties = FALSE) %>% ungroup()
top_combined <- bind_rows(top_up_combined, top_down_combined)

# Row order: UP terms first (grouped by cluster, desc NES), then DOWN terms
row_order_combined <- top_combined %>%
  arrange(desc(NES > 0), cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_combined <- expand_grid(
    pathway = unique(top_combined$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_combined)))

lim_combined <- max(abs(plot_df_combined$NES), na.rm = TRUE)

p_combined <- ggplot(plot_df_combined, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_combined, lim_combined),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "BD vs Control \n all enriched pathways", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

# ---- Save --------------------------------------------------------------------
write_csv(fgsea_up_df,   "fgsea_UP_df.csv")
write_csv(fgsea_down_df, "fgsea_DOWN_df.csv")

ggsave("tile_UP.pdf",       p_up,       width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_DOWN.pdf",     p_down,     width = 9,  height = 7,  device = cairo_pdf)
ggsave("tile_combined.pdf", p_combined, width = 10, height = 10, device = cairo_pdf)

# ---- Overall: top N per cluster by |NES|, ignoring direction -----------------
# fgsea_all_df already has results for all pathways in both directions.
# Just pick top N by |NES| per cluster regardless of sign.

top_overall <- fgsea_all_df %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall <- top_overall %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall <- expand_grid(
    pathway = unique(top_overall$pathway),
    cluster = unique(fgsea_all_df$cluster)
  ) %>%
  left_join(fgsea_all_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall)))

lim_overall <- max(abs(plot_df_overall$NES), na.rm = TRUE)

p_overall <- ggplot(plot_df_overall, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall, lim_overall),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall most affected pathways \nper cluster (top 5 by |NES|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

ggsave("tile_overall.pdf", p_overall, width = 10, height = 10, device = cairo_pdf)


# ---- FGSEA on merged UP+DOWN DE genes ranked by abs(t) -----------------------
# Question: what pathways are most dysregulated overall, ignoring direction?
# Rank by abs(t) — most dysregulated genes (either direction) at the top

all_DE <- bind_rows(Up_Regulated, Down_Regulated)

fgsea_overall_results <- lapply(clusters, function(cl) {
  ranked <- all_DE %>%
    filter(assay == cl) %>%
    arrange(desc(abs(t))) %>%
    { setNames(abs(.$t), .$ID) }

  fgsea(pathways    = all_gmt,
        stats       = ranked,
        minSize     = 5,
        maxSize     = 500,
        nPermSimple = NPERM) %>%
    as_tibble() %>%
    mutate(cluster = cl)
})

fgsea_overall_df <- bind_rows(fgsea_overall_results) %>%
  dplyr::select(-leadingEdge) %>%
  arrange(cluster, padj)

# Top N per cluster by NES
top_overall2 <- fgsea_overall_df %>%
  group_by(cluster) %>%
  slice_max(NES, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

row_order_overall2 <- top_overall2 %>%
  arrange(cluster, desc(NES)) %>%
  mutate(pathway = str_trunc(pathway, 50)) %>%
  pull(pathway) %>%
  unique()

plot_df_overall2 <- expand_grid(
    pathway = unique(top_overall2$pathway),
    cluster = unique(fgsea_overall_df$cluster)
  ) %>%
  left_join(fgsea_overall_df %>% dplyr::select(pathway, cluster, NES, padj),
            by = c("pathway", "cluster")) %>%
  mutate(pathway = str_trunc(pathway, 50),
         sig     = !is.na(padj) & padj < 0.05,
         pathway = factor(pathway, levels = rev(row_order_overall2)))

lim_overall2 <- max(abs(plot_df_overall2$NES), na.rm = TRUE)

p_overall2 <- ggplot(plot_df_overall2, aes(x = cluster, y = pathway, fill = NES)) +
  geom_tile(aes(color = sig), linewidth = 0.4, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim_overall2, lim_overall2),
                       na.value = "grey85", name = "NES") +
  coord_fixed() +
  labs(title = "Overall dysregulated pathways \n(ranked by |t|, UP+DOWN merged)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 7),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold"))

write_csv(fgsea_overall_df, "fgsea_overall_df.csv")
ggsave("tile_overall.pdf", p_overall2, width = 10, height = 8, device = cairo_pdf)
p_up
p_down
p_combined
p_overall
p_overall2

kegg_gmt <- parse_gmt("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.entrez.gmt")

library(dplyr)
library(purrr)
library(ggplot2)

# Pipeline-safe logging
log_msg <- function(msg, level = "INFO") {
  cat(sprintf("%s | %s | %s\n", Sys.time(), level, msg), file = stderr())
}

# 1. Parse GMT
parse_gmt <- function(gmt_path) {
  lines <- readLines(gmt_path)
  gmt_list <- strsplit(lines, "\t")
  names(gmt_list) <- sapply(gmt_list, `[[`, 1)
  lapply(gmt_list, function(x) x[c(-1, -2)]) 
}

# 2. Optimized ORA Function
run_manual_ora <- function(df, gmt_list, direction = "up", p_type = "adj.P.Val", p_cutoff = 0.05, top_n = 5) {
  universe <- unique(df$gene.name)
  
  # Safe subsetting to prevent dplyr context corruption
  if (direction == "up") {
    deg_df <- df %>% filter(logFC > 0 & !!sym(p_type) < p_cutoff)
  } else {
    deg_df <- df %>% filter(logFC < 0 & !!sym(p_type) < p_cutoff)
  }

  if (nrow(deg_df) == 0) {
    log_msg(sprintf("No genes passed thresholds for direction: %s", direction), "WARN")
    return(NULL)
  }

  # Bulletproof split bypassing the pipe
  split_dfs <- split(deg_df, deg_df[["assayID"]])
  
  # Run hypergeomtric per cluster
  results <- map_dfr(split_dfs, function(cluster_df) {
    cluster_genes <- unique(cluster_df$gene.name)
    k_total <- length(cluster_genes)
    curr_cluster <- cluster_df[["assayID"]][1]
    
    pathway_res <- imap_dfr(gmt_list, function(pathway_genes, pathway_name) {
      white_balls <- intersect(pathway_genes, universe)
      m <- length(white_balls)
      n <- length(universe) - m
      k <- length(intersect(cluster_genes, white_balls))
      
      p_val <- phyper(k - 1, m, n, k_total, lower.tail = FALSE)
      data.frame(term = pathway_name, p_val = p_val, overlap = k, pathway_size = m)
    })
    
    pathway_res %>%
      mutate(adj_p = p.adjust(p_val, method = "BH")) %>%
      filter(p_val < 0.05) %>% 
      arrange(p_val) %>%
      slice_head(n = top_n) %>%
      mutate(assayID = curr_cluster,
             mean_t = mean(cluster_df$t, na.rm = TRUE))
  })
  
  return(results)
}

# 3. Equal-Dimension Tile Plotter
plot_enrichment_tiles <- function(res_df, title) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)

  ggplot(res_df, aes(x = assayID, y = term, fill = mean_t)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    coord_fixed() + # Forces identical width/height for tiles
    theme_minimal() +
    labs(title = title, x = "Cell Cluster", y = "Enriched Process", fill = "Avg t-stat") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# --- Execution ---
gmt_path <- "../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt"
gmt_gobp <- parse_gmt(gmt_path)

up_enrich_df   <- run_manual_ora(AllDF, gmt_gobp, direction = "up", p_type = "P.Value")
down_enrich_df <- run_manual_ora(AllDF, gmt_gobp, direction = "down", p_type = "P.Value")

p_up <- plot_enrichment_tiles(up_enrich_df, "BD Case Enrichment (Up)")
p_down <- plot_enrichment_tiles(down_enrich_df, "Control Enrichment (Down)")

up_enrich_df

plot_enrichment_tiles <- function(res_df, title) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)

  ggplot(res_df, aes(x = assayID, y = term, fill = mean_t)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    coord_fixed() + # Forces perfect squares
    theme_minimal() +
    labs(title = title, x = "Cell Cluster", y = "Enriched Process", fill = "Avg t-stat") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))
}

p_up <- plot_enrichment_tiles(up_enrich_df, "BD Case Enrichment (Up)")
p_down <- plot_enrichment_tiles(down_enrich_df, "Control Enrichment (Down)")

p_up

# ─────────────────────────────────────────────
# RUN — adjust toggles here
# ─────────────────────────────────────────────

CHOOSE_METHOD    <- "degs2_only"   # "union", "intersect", "degs1_only", "degs2_only"
CHOSEN_STAT      <- "t"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 10
USE_ADJ_ORG      <- FALSE
USE_ADJ_TISS     <- TRUE
PVAL_THRES_ORG   <- 0.05
PVAL_THRES_TISS  <- 0.05
LFC_THRESH       <- 0
SAVE_PLOTS       <- FALSE
FIGURE_PATH      <- "./figures/"
PLOT_WIDTH       <- 7
PLOT_HEIGHT      <- 5
PLOT_RES         <- 200

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

# --- 1. Feature selection ---
final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = CHOOSE_METHOD,
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  use_sig_filter  = TRUE,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

# --- 2. Matrix prep (t-statistic) ---
mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop=FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop=FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse="x"),
                 paste(dim(mat_tiss_final), collapse="x")))

# --- 3. Spearman correlation ---
cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)
correlation_qc_report(cor_result)


plot_subtitle <- sprintf(
  "Organoid: nominal p<%.2f | Tissue: nominal p<%.2f \n Genes: %d",
  PVAL_THRES_ORG, PVAL_THRES_TISS, length(final_features)
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = "Organoid vs Tissue — Spearman (t-statistic)",
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  save        = SAVE_PLOTS,
  file_prefix = "BD_spearman",
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# --- 4. Scatter plots ---
library(ggplot2)
library(ggrepel)

# Single pair — no gene labels
plot_scatter_concordance(
  mat_org      = mat_org_final,
  mat_tiss     = mat_tiss_final,
  org_cluster  = "Astro",
  tiss_cluster = "Astro",
  stat_label   = "t-statistic",
  org_label    = "Organoid",
  tiss_label   = "Tissue"
)

# All biologically expected pairs — t-statistic
plot_scatter_all_pairs(
  mat_org     = mat_org_final,
  mat_tiss    = mat_tiss_final,
  pairs       = list(
    c("Astro",      "Astro"),
    c("Neuron",     "IN"),
    c("Neuron",     "IN"),
    c("IPC",        "OPC"),
    c("Fibroblast", "Mural")
  ),
  stat_label  = "t-statistic",
  org_label   = "Organoid",
  tiss_label  = "Tissue",
  save        = SAVE_PLOTS,
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# --- 5. logFC scatter (separate matrix prep) ---
mat_org_lfc  <- prep_cluster_matrix(org_df_flat,  final_features, "logFC")
mat_tiss_lfc <- prep_cluster_matrix(tiss_df_flat, final_features, "logFC")

common_genes_lfc  <- sort(intersect(rownames(mat_org_lfc), rownames(mat_tiss_lfc)))
mat_org_lfc_final  <- as.matrix(mat_org_lfc[common_genes_lfc,  , drop=FALSE])
mat_tiss_lfc_final <- as.matrix(mat_tiss_lfc[common_genes_lfc, , drop=FALSE])

plot_scatter_all_pairs(
  mat_org     = mat_org_lfc_final,
  mat_tiss    = mat_tiss_lfc_final,
  pairs       = list(
    c("Astro",  "Astro"),
    c("Neuron", "IN"),
    c("Neuron", "EN")
  ),
  stat_label  = "logFC",
  org_label   = "Organoid",
  tiss_label  = "Tissue",
  save        = SAVE_PLOTS,
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# --- 6. Scatter with concordant gene labels (run after get_concordant_genes) ---
# Uncomment once concordant_results is available:
# conc_genes_astro <- concordant_results[["Astro_Astro"]]$gene
# plot_scatter_concordance(
#   mat_org         = mat_org_final,
#   mat_tiss        = mat_tiss_final,
#   org_cluster     = "Astro",
#   tiss_cluster    = "Astro",
#   highlight_genes = head(conc_genes_astro, 20),
#   stat_label      = "t-statistic",
#   save            = SAVE_PLOTS,
#   figure_path     = FIGURE_PATH
# )

# ─────────────────────────────────────────────
# RUN — adjust toggles here
# ─────────────────────────────────────────────

# Toggles
CHOSEN_STAT      <- "logFC"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 2
USE_ADJ_ORG      <- FALSE     # organoid: nominal p-value
USE_ADJ_TISS     <- FALSE     # tissue: nominal p-value (BD effect too small for FDR)
PVAL_THRES_ORG   <- 0.05
PVAL_THRES_TISS  <- 0.05
LFC_THRESH       <- 0         # raise to 0.25/0.5 for stricter DEG definition
SAVE_PLOTS       <- FALSE     # set TRUE to save all figures
FIGURE_PATH      <- "./figures/"
PLOT_WIDTH       <- 7
PLOT_HEIGHT      <- 5
PLOT_RES         <- 200

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

# --- 1. Feature selection ---
final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = "union",
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  use_sig_filter  = TRUE,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

# --- 2. Matrix prep ---
mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop=FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop=FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse="x"),
                 paste(dim(mat_tiss_final), collapse="x")))

# --- 3. Spearman correlation ---
cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)

correlation_qc_report(cor_result)
plot_subtitle <- sprintf(
  "Organoid: nominal p<%.2f | Tissue: nominal p<%.2f \n Genes: %d",
  PVAL_THRES_ORG, PVAL_THRES_TISS, length(final_features)
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = paste("Organoid vs Tissue — Spearman (",CHOSEN_STAT,")"),
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  save        = SAVE_PLOTS,
  file_prefix = "BD_spearman",
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# ─────────────────────────────────────────────
# RUN — adjust toggles here
# ─────────────────────────────────────────────

# Toggles
CHOOSE_METHOD    <- "intersect"
CHOSEN_STAT      <- "logFC"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 2
USE_ADJ_ORG      <- FALSE     # organoid: nominal p-value
USE_ADJ_TISS     <- FALSE     # tissue: nominal p-value (BD effect too small for FDR)
PVAL_THRES_ORG   <- 0.05
PVAL_THRES_TISS  <- 0.05
LFC_THRESH       <- 0         # raise to 0.25/0.5 for stricter DEG definition
SAVE_PLOTS       <- FALSE     # set TRUE to save all figures
FIGURE_PATH      <- "./figures/"
PLOT_WIDTH       <- 7
PLOT_HEIGHT      <- 5
PLOT_RES         <- 200

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

# --- 1. Feature selection ---
final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = CHOOSE_METHOD,
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  use_sig_filter  = TRUE,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

# --- 2. Matrix prep ---
mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop=FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop=FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse="x"),
                 paste(dim(mat_tiss_final), collapse="x")))

# --- 3. Spearman correlation ---
cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)

correlation_qc_report(cor_result)
plot_subtitle <- sprintf(
  "Organoid: nominal p<%.2f | Tissue: nominal p<%.2f \n Genes: %d",
  PVAL_THRES_ORG, PVAL_THRES_TISS, length(final_features)
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = paste("Organoid vs Tissue — Spearman (",CHOSEN_STAT,")"),
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  save        = SAVE_PLOTS,
  file_prefix = "BD_spearman",
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# ─────────────────────────────────────────────
# RUN — adjust toggles here
# ─────────────────────────────────────────────

# Toggles
CHOOSE_METHOD    <- "union"
CHOSEN_STAT      <- "logFC"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 2
USE_ADJ_ORG      <- FALSE     # organoid: nominal p-value
USE_ADJ_TISS     <- FALSE     # tissue: nominal p-value (BD effect too small for FDR)
PVAL_THRES_ORG   <- 0.05
PVAL_THRES_TISS  <- 0.05
LFC_THRESH       <- 0         # raise to 0.25/0.5 for stricter DEG definition
SAVE_PLOTS       <- FALSE     # set TRUE to save all figures
FIGURE_PATH      <- "./figures/"
PLOT_WIDTH       <- 7
PLOT_HEIGHT      <- 5
PLOT_RES         <- 200

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

# --- 1. Feature selection ---
final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = CHOOSE_METHOD,
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  use_sig_filter  = TRUE,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

# --- 2. Matrix prep ---
mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop=FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop=FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse="x"),
                 paste(dim(mat_tiss_final), collapse="x")))

# --- 3. Spearman correlation ---
cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)

correlation_qc_report(cor_result)
plot_subtitle <- sprintf(
  "Organoid: nominal p<%.2f | Tissue: nominal p<%.2f \n Genes: %d",
  PVAL_THRES_ORG, PVAL_THRES_TISS, length(final_features)
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = paste("Organoid vs Tissue — Spearman (",CHOSEN_STAT,")"),
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  save        = SAVE_PLOTS,
  file_prefix = "BD_spearman",
  figure_path = FIGURE_PATH,
  plot.width  = PLOT_WIDTH,
  plot.height = PLOT_HEIGHT,
  plot.res    = PLOT_RES
)

# Diagnose what's actually in the matrices
cat("=== Feature set size ===\n")
cat(sprintf("Final features: %d genes\n", length(final_features)))

cat("\n=== Matrix dimensions ===\n")
cat(sprintf("Organoid matrix:  %d genes x %d clusters\n", 
            nrow(mat_org_final), ncol(mat_org_final)))
cat(sprintf("Tissue matrix:    %d genes x %d clusters\n", 
            nrow(mat_tiss_final), ncol(mat_tiss_final)))

cat("\n=== Organoid matrix sample (first 5 genes x all clusters) ===\n")
print(round(mat_org_final[1:5, ], 3))

cat("\n=== Tissue matrix sample (first 5 genes x all clusters) ===\n")
print(round(mat_tiss_final[1:5, ], 3))

cat("\n=== NA counts per cluster ===\n")
cat("Organoid NAs per cluster:\n")
print(colSums(is.na(mat_org_final)))
cat("Tissue NAs per cluster:\n")
print(colSums(is.na(mat_tiss_final)))

cat("\n=== t-stat range check ===\n")
cat("Organoid t range:", round(range(mat_org_final, na.rm=TRUE), 3), "\n")
cat("Tissue t range:  ", round(range(mat_tiss_final, na.rm=TRUE), 3), "\n")

cat("\n=== Spot check: Astro vs Astro correlation ===\n")
org_astro  <- mat_org_final[, "Astro"]
tiss_astro <- mat_tiss_final[, "Astro"]
keep <- !is.na(org_astro) & !is.na(tiss_astro)
cat(sprintf("Shared non-NA genes: %d\n", sum(keep)))
cat(sprintf("Spearman r: %.4f\n", 
            cor(org_astro[keep], tiss_astro[keep], method="spearman")))

# Check if the issue is sign flipping — maybe BD vs CTR is coded differently
cat("=== Direction check: mean t per cluster ===\n")
cat("Organoid (positive = up in BD?):\n")
print(round(colMeans(mat_org_final, na.rm=TRUE), 4))

cat("\nTissue (positive = up in BD?):\n")
print(round(colMeans(mat_tiss_final, na.rm=TRUE), 4))

# Check contrast coding in original dreamlet results
cat("\n=== Organoid contrast check ===\n")
cat("coefNames or contrast used:\n")
print(unique(org_df_flat$contrast))   # adjust column name if different

cat("\n=== Tissue contrast check ===\n")
print(unique(tiss_df_flat$contrast))  # adjust column name if different

# Visual check — scatter plot Astro vs Astro t-stats
org_astro  <- mat_org_final[, "Astro"]
tiss_astro <- mat_tiss_final[, "Astro"]
keep       <- !is.na(org_astro) & !is.na(tiss_astro)

plot(
  org_astro[keep], tiss_astro[keep],
  xlab = "Organoid Astro t-stat (BD vs CTR)",
  ylab = "Tissue Astro t-stat (BD vs CTR)",
  main = sprintf("Astro BD signature concordance\nr=%.3f, n=%d genes",
                 cor(org_astro[keep], tiss_astro[keep], method="spearman"),
                 sum(keep)),
  pch  = 16, cex = 0.4, col = "#00000040"
)
abline(h=0, v=0, col="red", lty=2)
abline(lm(tiss_astro[keep] ~ org_astro[keep]), col="blue", lwd=1.5)



# Quick comparison: t vs z.std correlation for Astro vs Astro
for (stat in c("t", "z.std", "logFC")) {
  mat_o <- prep_cluster_matrix(org_df_flat,  final_features, stat)
  mat_ti <- prep_cluster_matrix(tiss_df_flat, final_features, stat)
  cg <- sort(intersect(rownames(mat_o), rownames(mat_ti)))
  x  <- as.matrix(mat_o[cg,  , drop=FALSE])[, "Astro"]
  y  <- as.matrix(mat_ti[cg, , drop=FALSE])[, "Astro"]
  keep <- !is.na(x) & !is.na(y)
  cat(sprintf("%-8s — Astro vs Astro r=%.4f (n=%d)\n",
              stat,
              cor(x[keep], y[keep], method="spearman"),
              sum(keep)))
}

# Run cor.test for ALL biologically expected pairs and report properly
expected_pairs <- list(
  c("Astro",      "Astro"),
  c("Neuron",     "IN"),
  c("Neuron",     "EN"),
  c("IPC",        "OPC"),
  c("Fibroblast", "Mural"),
  c("Fibroblast", "Endo")
)

cat("=== Correlation significance for expected pairs ===\n")
cat(sprintf("%-12s %-10s %8s %12s %6s\n",
            "Org", "Tissue", "r", "p-value", "n"))
cat(strrep("-", 55), "\n")

results_df <- map_dfr(expected_pairs, function(pair) {
  org_ct  <- pair[1]
  tiss_ct <- pair[2]

  x    <- mat_org_final[,  org_ct]
  y    <- mat_tiss_final[, tiss_ct]
  keep <- !is.na(x) & !is.na(y)
  n    <- sum(keep)

  if (n < 10) {
    cat(sprintf("%-12s %-10s %8s %12s %6d  [MASKED]\n",
                org_ct, tiss_ct, "NA", "NA", n))
    return(NULL)
  }

  test <- cor.test(x[keep], y[keep], method = "spearman")
  r    <- round(test$estimate, 4)
  pval <- test$p.value

  cat(sprintf("%-12s %-10s %8.4f %12.3e %6d\n",
              org_ct, tiss_ct, r, pval, n))

  data.frame(
    organoid   = org_ct,
    tissue     = tiss_ct,
    r          = r,
    pval       = pval,
    n_genes    = n
  )
})

# FDR across expected pairs
results_df$padj <- p.adjust(results_df$pval, method = "BH")

cat("\n=== With BH correction across expected pairs ===\n")
print(results_df %>%
  mutate(
    sig    = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  ) %>%
  arrange(padj),
  row.names = FALSE
)


get_concordant_genes <- function(
    mat_org, mat_tiss,
    org_df, tiss_df,
    org_cluster,
    tiss_cluster,
    pval_thres_org  = 0.05,
    pval_thres_tiss = 0.05,
    top_n           = Inf
) {
  x    <- mat_org[,  org_cluster]
  y    <- mat_tiss[, tiss_cluster]
  keep <- !is.na(x) & !is.na(y)
  genes <- names(x)[keep]

  # Significant in organoid — nominal
  sig_org <- org_df %>%
    filter(assay == org_cluster, P.Value < pval_thres_org) %>%
    pull(ID) %>%
    unique()

  # Significant in tissue — nominal
  sig_tiss <- tiss_df %>%
    filter(assay == tiss_cluster, P.Value < pval_thres_tiss) %>%
    pull(ID) %>%
    unique()

  # Concordant: significant in both + same direction
  concordant <- genes[
    genes %in% sig_org  &
    genes %in% sig_tiss &
    sign(x[genes]) == sign(y[genes])
  ]

  if (length(concordant) == 0) {
    cat("  No concordant genes found.\n")
    return(data.frame())
  }

  result <- data.frame(
    gene      = concordant,
    t_org     = round(x[concordant], 3),
    t_tiss    = round(y[concordant], 3),
    mean_abs  = round((abs(x[concordant]) + abs(y[concordant])) / 2, 3),
    direction = ifelse(x[concordant] > 0, "UP_BD", "DOWN_BD")
  ) %>%
    arrange(desc(mean_abs))

  if (is.finite(top_n)) result <- head(result, top_n)

  return(result)
}

# Full concordant gene analysis across all significant pairs
concordant_results <- list()

pairs_to_extract <- list(
  list(org="Astro",      tiss="Astro", label="Astro_Astro"),
  list(org="Neuron",     tiss="IN",    label="Neuron_IN"),
  list(org="Fibroblast", tiss="Mural", label="Fibroblast_Mural")
)

for (pair in pairs_to_extract) {
  cat(sprintf("\n=== Concordant genes: %s ===\n", pair$label))

  conc <- get_concordant_genes(
    mat_org_final, mat_tiss_final,
    org_df_flat,   tiss_df_flat,
    org_cluster    = pair$org,
    tiss_cluster   = pair$tiss,
    top_n          = Inf   # return all
  )

  cat(sprintf("Total concordant:  %d\n", nrow(conc)))
  cat(sprintf("UP in BD:          %d\n", sum(conc$direction == "UP_BD")))
  cat(sprintf("DOWN in BD:        %d\n", sum(conc$direction == "DOWN_BD")))
  cat("Top 15:\n")
  print(head(conc, 15), row.names=FALSE)

  concordant_results[[pair$label]] <- conc
}

# Save all concordant gene lists
saveRDS(concordant_results, "BD_concordant_genes_organoid_tissue.rds")

# Summary table across pairs
cat("\n=== SUMMARY: Concordant gene counts ===\n")
summary_df <- map_dfr(names(concordant_results), function(label) {
  df <- concordant_results[[label]]
  data.frame(
    pair       = label,
    n_total    = nrow(df),
    n_up_BD    = sum(df$direction == "UP_BD"),
    n_down_BD  = sum(df$direction == "DOWN_BD"),
    top_gene   = df$gene[1]
  )
})
print(summary_df, row.names=FALSE)


# Save concordant gene lists for enrichment
write.csv(concordant_results$Astro_Astro,
          "concordant_BD_Astro_organoid_tissue.csv",
          row.names=FALSE)
write.csv(concordant_results$Neuron_IN,
          "concordant_BD_Neuron_IN_organoid_tissue.csv",
          row.names=FALSE)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(purrr)

# ─────────────────────────────────────────────
# 1. UTILITIES
# ─────────────────────────────────────────────

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
  flush(stderr())
}

# ─────────────────────────────────────────────
# 2. FEATURE SELECTION
# ─────────────────────────────────────────────

get_feature_genes <- function(
    df_org, df_tiss,
    method          = "union",
    stat_type       = "t",
    use_concordance = FALSE,
    pval_thres_org  = 0.05,
    pval_thres_tiss = 0.05,
    use_adj_org     = FALSE,
    use_adj_tiss    = TRUE,
    min_n_pairs     = 10
) {
  # 1. Background: genes present in both datasets
  common_background <- intersect(unique(df_org$ID), unique(df_tiss$ID))
  log_info(sprintf("Background genes (shared): %d", length(common_background)))

  # 2. Organoid significance — nominal by default (underpowered for FDR)
  if (use_adj_org) {
    if (!"adj.P.Val" %in% colnames(df_org)) stop("adj.P.Val not found in df_org")
    sig_ids_org <- unique(df_org$ID[df_org$adj.P.Val < pval_thres_org])
    log_info(sprintf("Organoid: %d genes at adj.P.Val < %.2f",
                     length(sig_ids_org), pval_thres_org))
  } else {
    sig_ids_org <- unique(df_org$ID[df_org$P.Value < pval_thres_org])
    log_info(sprintf("Organoid: %d genes at nominal P.Value < %.2f",
                     length(sig_ids_org), pval_thres_org))
  }

  # 3. Tissue significance — FDR by default (well-powered)
  if (use_adj_tiss) {
    if (!"adj.P.Val" %in% colnames(df_tiss)) stop("adj.P.Val not found in df_tiss")
    sig_ids_tiss <- unique(df_tiss$ID[df_tiss$adj.P.Val < pval_thres_tiss])
    log_info(sprintf("Tissue: %d genes at adj.P.Val < %.2f",
                     length(sig_ids_tiss), pval_thres_tiss))
  } else {
    sig_ids_tiss <- unique(df_tiss$ID[df_tiss$P.Value < pval_thres_tiss])
    log_info(sprintf("Tissue: %d genes at nominal P.Value < %.2f",
                     length(sig_ids_tiss), pval_thres_tiss))
  }

  # 4. Union or intersection
  if (method == "union") {
    target_sig <- union(sig_ids_org, sig_ids_tiss)
  } else {
    target_sig <- intersect(sig_ids_org, sig_ids_tiss)
  }

  candidate_genes <- intersect(target_sig, common_background)
  log_info(sprintf("Candidate genes after %s + background filter: %d",
                   method, length(candidate_genes)))

  # 5. Optional concordance — per cell type majority vote
  # Avoids averaging t/logFC across cell types which loses specificity
  if (use_concordance) {
    trend_org <- df_org %>%
      filter(ID %in% candidate_genes) %>%
      group_by(ID, assay) %>%
      summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
      rename(val_org = val)

    trend_tiss <- df_tiss %>%
      filter(ID %in% candidate_genes) %>%
      group_by(ID, assay) %>%
      summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
      rename(val_tiss = val)

    concordance <- inner_join(trend_org, trend_tiss, by = c("ID", "assay")) %>%
      mutate(concordant = sign(val_org) == sign(val_tiss)) %>%
      group_by(ID) %>%
      summarize(
        n_types      = n(),
        n_concordant = sum(concordant),
        pct_conc     = n_concordant / n_types,
        .groups      = "drop"
      ) %>%
      filter(pct_conc >= 0.5)

    candidate_genes <- concordance$ID
    log_info(sprintf("After majority concordance filter: %d genes",
                     length(candidate_genes)))
  }

  return(candidate_genes)
}

# ─────────────────────────────────────────────
# 3. MATRIX PREPARATION
# ─────────────────────────────────────────────

prep_cluster_matrix <- function(df, genes, stat_type) {
  df %>%
    as.data.frame() %>%
    filter(ID %in% genes) %>%
    group_by(assay, ID) %>%
    summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = assay, values_from = val) %>%
    column_to_rownames("ID")
}

# ─────────────────────────────────────────────
# 4. CORRELATION WITH SIGNIFICANCE + N MASKING
# ─────────────────────────────────────────────

get_pairwise_n_matrix <- function(mat_x, mat_y) {
  indicator_x <- !is.na(mat_x)
  indicator_y <- !is.na(mat_y)
  as.matrix(t(indicator_x) %*% indicator_y)
}

compute_correlation_with_pval <- function(
    mat_org, mat_tiss,
    method = "spearman",
    min_n  = 10
) {
  n_org  <- ncol(mat_org)
  n_tiss <- ncol(mat_tiss)

  cor_mat  <- matrix(NA, nrow = n_org, ncol = n_tiss,
                     dimnames = list(colnames(mat_org), colnames(mat_tiss)))
  pval_mat <- matrix(NA, nrow = n_org, ncol = n_tiss,
                     dimnames = list(colnames(mat_org), colnames(mat_tiss)))
  n_mat    <- get_pairwise_n_matrix(mat_org, mat_tiss)

  for (i in seq_len(n_org)) {
    for (j in seq_len(n_tiss)) {
      x    <- mat_org[, i]
      y    <- mat_tiss[, j]
      keep <- !is.na(x) & !is.na(y)
      n    <- sum(keep)

      if (n < min_n) {
        log_info(sprintf(
          "Masking pair [%s x %s]: only %d shared genes (min_n=%d)",
          colnames(mat_org)[i], colnames(mat_tiss)[j], n, min_n
        ))
        next
      }

      test_result   <- cor.test(x[keep], y[keep], method = method)
      cor_mat[i, j] <- test_result$estimate
      pval_mat[i, j] <- test_result$p.value
    }
  }

  # FDR correction across all pairs
  pval_adj_mat <- matrix(
    p.adjust(as.vector(pval_mat), method = "BH"),
    nrow = n_org, ncol = n_tiss,
    dimnames = list(colnames(mat_org), colnames(mat_tiss))
  )

  list(
    cor  = cor_mat,
    pval = pval_mat,
    padj = pval_adj_mat,
    n    = n_mat
  )
}

# ─────────────────────────────────────────────
# 5. VISUALIZATION
# ─────────────────────────────────────────────
plot_correlation_heatmap <- function(
    cor_mat,
    pval_mat,
    n_mat,
    title       = "Organoid vs Tissue",
    subtitle    = "",
    sig_only    = FALSE,
    padj_thresh = 0.05
) {
  # Replace NA with 0 ONLY for clustering
  mat_for_clustering <- cor_mat
  mat_for_clustering[is.na(mat_for_clustering)] <- 0

  display_mat <- cor_mat
  if (sig_only) {
    display_mat[pval_mat > padj_thresh] <- NA
  }

  # Significance stars
  sig_stars <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                      dimnames = dimnames(pval_mat))
  sig_stars[!is.na(pval_mat) & pval_mat < 0.05]  <- "*"
  sig_stars[!is.na(pval_mat) & pval_mat < 0.01]  <- "**"
  sig_stars[!is.na(pval_mat) & pval_mat < 0.001] <- "***"

  # Cell annotation
  cell_labels <- matrix("NA", nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                        dimnames = dimnames(cor_mat))
  not_na <- !is.na(cor_mat)
  cell_labels[not_na] <- sprintf("%.2f%s", cor_mat[not_na], sig_stars[not_na])

  limit <- max(abs(cor_mat), na.rm = TRUE)

  # Pre-compute clustering on NA-imputed matrix
  row_order <- hclust(dist(mat_for_clustering))$order
  col_order <- hclust(dist(t(mat_for_clustering)))$order

  # Combine title + subtitle into main — pheatmap supports \n
  full_title <- if (nchar(subtitle) > 0) paste0(title, "\n", subtitle) else title

  pheatmap(
    display_mat[row_order, col_order],   # reorder manually
    display_numbers  = cell_labels[row_order, col_order],
    number_color     = "black",
    color            = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
    breaks           = seq(-limit, limit, length.out = 101),
    border_color     = "grey90",
    na_col           = "grey85",
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    fontsize         = 11,
    fontsize_number  = 8,
    main             = full_title,
    legend_breaks    = c(-limit, -limit/2, 0, limit/2, limit),
    legend_labels    = c(
      sprintf("%.2f", -limit),
      sprintf("%.2f", -limit/2),
      "0",
      sprintf("%.2f",  limit/2),
      sprintf("%.2f",  limit)
    ),
    angle_col = 45
  )
}

# ─────────────────────────────────────────────
# 6. QC REPORT
# ─────────────────────────────────────────────

correlation_qc_report <- function(cor_result) {
  cat("\n", strrep("=", 60), "\n")
  cat("  CORRELATION QC REPORT\n")
  cat(strrep("=", 60), "\n\n")

  cor_vals <- as.vector(cor_result$cor)
  cor_vals <- cor_vals[!is.na(cor_vals)]

  cat(sprintf("  Pairs computed:        %d of %d\n",
              sum(!is.na(cor_result$cor)), length(cor_result$cor)))
  cat(sprintf("  Pairs masked (low N):  %d\n",
              sum(is.na(cor_result$cor))))
  cat(sprintf("  Mean r:                %.3f\n", mean(cor_vals)))
  cat(sprintf("  Std r:                 %.3f\n", sd(cor_vals)))
  cat(sprintf("  Min r:                 %.3f\n", min(cor_vals)))
  cat(sprintf("  Max r:                 %.3f\n", max(cor_vals)))
  cat(sprintf("  Sig pairs (padj<0.05): %d of %d (%.1f%%)\n",
              sum(cor_result$padj < 0.05, na.rm = TRUE),
              sum(!is.na(cor_result$padj)),
              100 * mean(cor_result$padj < 0.05, na.rm = TRUE)))

  cat("\n[Best Match per Organoid Cluster — Spearman]\n")
  best_df <- data.frame(
    organoid    = rownames(cor_result$cor),
    best_tissue = colnames(cor_result$cor)[
      apply(cor_result$cor, 1, function(x) which.max(x))
    ],
    max_r = round(apply(cor_result$cor, 1, max, na.rm = TRUE), 3),
    min_r = round(apply(cor_result$cor, 1, min, na.rm = TRUE), 3)
  ) %>%
    mutate(delta = round(max_r - min_r, 3)) %>%
    arrange(desc(max_r))

  print(best_df, row.names = FALSE)
  cat("\n  [delta = discriminability: higher = cleaner mapping]\n")
  cat(strrep("=", 60), "\n\n")
}

COEF <- "Diff"
organoid_df <- topTable(readRDS('2026_03_06_Scanpy_020_2dpass_annotated_CLASS_cova_DiffPDvsCTRL_1.9.rds'), coef = COEF, number = Inf)
tissue_df <- topTable(readRDS('Lyra_data/dreamlet_class_upitt.rds'), coef = COEF, number = Inf)
#tissue_df <- readRDS('Lyra_data/subclass_upitt_dreamlet_toptable.rds')
tissue_sets <- get_gene_list(tissue_df, "T_")
organoid_sets <- get_gene_list(organoid_df, "O_")

# ─────────────────────────────────────────────
# 7. RUN
# ─────────────────────────────────────────────

# Toggles
CHOSEN_METHOD    <- "intersect"
CHOSEN_STAT      <- "t"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 10
USE_ADJ_ORG      <- FALSE   # nominal — organoid underpowered for FDR
USE_ADJ_TISS     <- FALSE    # FDR   — tissue well-powered
PVAL_THRES_ORG   <- 0.01
PVAL_THRES_TISS  <- 0.01

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

# Feature selection
final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = CHOSEN_METHOD,
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

# Matrix prep
mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse = "x"),
                 paste(dim(mat_tiss_final), collapse = "x")))

# Correlation
cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)

# QC report
correlation_qc_report(cor_result)

# Plot
plot_subtitle <- sprintf(
  "Method: %s | Stat: %s | Concordance: %s | Genes: %d\nOrganoid: %s p<%.2f | Tissue: %s p<%.2f",
  toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),
  ifelse(USE_ADJ_ORG,  "FDR", "nominal"), PVAL_THRES_ORG,
  ifelse(USE_ADJ_TISS, "FDR", "nominal"), PVAL_THRES_TISS
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = "Organoid vs Tissue — Spearman Correlation",
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  padj_thresh = 0.05
)

cat("=== Organoid t-stat summary per assay ===\n")
org_df_flat %>%
  group_by(assay) %>%
  summarize(
    n      = n(),
    mean_t = round(mean(t, na.rm=TRUE), 3),
    sd_t   = round(sd(t, na.rm=TRUE), 3),
    min_t  = round(min(t, na.rm=TRUE), 3),
    max_t  = round(max(t, na.rm=TRUE), 3),
    n_sig  = sum(P.Value < 0.05, na.rm=TRUE),
    .groups = "drop"
  ) %>% print()

cat("\n=== Tissue t-stat summary per assay ===\n")
tiss_df_flat %>%
  group_by(assay) %>%
  summarize(
    n      = n(),
    mean_t = round(mean(t, na.rm=TRUE), 3),
    sd_t   = round(sd(t, na.rm=TRUE), 3),
    min_t  = round(min(t, na.rm=TRUE), 3),
    max_t  = round(max(t, na.rm=TRUE), 3),
    n_sig  = sum(adj.P.Val < 0.05, na.rm=TRUE),
    .groups = "drop"
  ) %>% print()

cat("\n=== Shared significant genes per matching pair ===\n")
matching_pairs <- list(
  Astro      = "Astro",
  Neuron     = "IN",
  IPC        = "OPC",
  Fibroblast = "Mural"
)
for (org_ct in names(matching_pairs)) {
  tiss_ct    <- matching_pairs[[org_ct]]
  org_genes  <- org_df_flat$ID[
    org_df_flat$assay == org_ct & org_df_flat$P.Value < 0.05
  ]
  tiss_genes <- tiss_df_flat$ID[
    tiss_df_flat$assay == tiss_ct & tiss_df_flat$adj.P.Val < 0.05
  ]
  overlap <- intersect(org_genes, tiss_genes)
  cat(sprintf("  %s <-> %s | org: %d | tiss: %d | overlap: %d\n",
              org_ct, tiss_ct,
              length(org_genes), length(tiss_genes), length(overlap)))
}

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(purrr)
library(gridExtra)

# ─────────────────────────────────────────────
# 1. UTILITIES
# ─────────────────────────────────────────────

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
  flush(stderr())
}

# ─────────────────────────────────────────────
# 2. FEATURE SELECTION
# ─────────────────────────────────────────────

get_feature_genes <- function(
    df_org, df_tiss,
    method          = "union",
    stat_type       = "t",
    use_concordance = FALSE,
    pval_thres_org  = 0.05,
    pval_thres_tiss = 0.05,
    use_adj_org     = FALSE,
    use_adj_tiss    = FALSE,
    use_sig_filter  = TRUE,
    min_n_pairs     = 10
) {
  # 1. Background: genes present in both datasets
  common_background <- intersect(unique(df_org$ID), unique(df_tiss$ID))
  log_info(sprintf("Background genes (shared): %d", length(common_background)))

  if (!use_sig_filter) {
    log_info("Significance filter DISABLED — using all shared background genes")
    candidate_genes <- common_background

  } else {
    # Organoid significance
    if (use_adj_org) {
      if (!"adj.P.Val" %in% colnames(df_org)) stop("adj.P.Val not found in df_org")
      sig_ids_org <- unique(df_org$ID[df_org$adj.P.Val < pval_thres_org])
      log_info(sprintf("Organoid: %d genes at adj.P.Val < %.2f",
                       length(sig_ids_org), pval_thres_org))
    } else {
      sig_ids_org <- unique(df_org$ID[df_org$P.Value < pval_thres_org])
      log_info(sprintf("Organoid: %d genes at nominal P.Value < %.2f",
                       length(sig_ids_org), pval_thres_org))
    }

    # Tissue significance
    if (use_adj_tiss) {
      if (!"adj.P.Val" %in% colnames(df_tiss)) stop("adj.P.Val not found in df_tiss")
      sig_ids_tiss <- unique(df_tiss$ID[df_tiss$adj.P.Val < pval_thres_tiss])
      log_info(sprintf("Tissue: %d genes at adj.P.Val < %.2f",
                       length(sig_ids_tiss), pval_thres_tiss))
    } else {
      sig_ids_tiss <- unique(df_tiss$ID[df_tiss$P.Value < pval_thres_tiss])
      log_info(sprintf("Tissue: %d genes at nominal P.Value < %.2f",
                       length(sig_ids_tiss), pval_thres_tiss))
    }

    # Union or intersection of significant genes
    if (method == "union") {
      target_sig <- union(sig_ids_org, sig_ids_tiss)
    } else {
      target_sig <- intersect(sig_ids_org, sig_ids_tiss)
    }

    candidate_genes <- intersect(target_sig, common_background)
    log_info(sprintf("Candidate genes after %s + background: %d",
                     method, length(candidate_genes)))
  }

  # Concordance — per cell type majority vote
  if (use_concordance) {
    trend_org <- df_org %>%
      filter(ID %in% candidate_genes) %>%
      group_by(ID, assay) %>%
      summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
      rename(val_org = val)

    trend_tiss <- df_tiss %>%
      filter(ID %in% candidate_genes) %>%
      group_by(ID, assay) %>%
      summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
      rename(val_tiss = val)

    concordance <- inner_join(trend_org, trend_tiss, by = c("ID", "assay")) %>%
      mutate(concordant = sign(val_org) == sign(val_tiss)) %>%
      group_by(ID) %>%
      summarize(
        n_types      = n(),
        n_concordant = sum(concordant),
        pct_conc     = n_concordant / n_types,
        .groups      = "drop"
      ) %>%
      filter(pct_conc >= 0.5)

    candidate_genes <- concordance$ID
    log_info(sprintf("After majority concordance filter: %d genes",
                     length(candidate_genes)))
  }

  return(candidate_genes)
}

# ─────────────────────────────────────────────
# 3. MATRIX PREPARATION
# ─────────────────────────────────────────────

prep_cluster_matrix <- function(df, genes, stat_type) {
  df %>%
    as.data.frame() %>%
    filter(ID %in% genes) %>%
    group_by(assay, ID) %>%
    summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = assay, values_from = val) %>%
    column_to_rownames("ID")
}

# ─────────────────────────────────────────────
# 4. CORRELATION WITH SIGNIFICANCE + N MASKING
# ─────────────────────────────────────────────

get_pairwise_n_matrix <- function(mat_x, mat_y) {
  indicator_x <- !is.na(mat_x)
  indicator_y <- !is.na(mat_y)
  as.matrix(t(indicator_x) %*% indicator_y)
}

compute_correlation_with_pval <- function(
    mat_org, mat_tiss,
    method = "spearman",
    min_n  = 10
) {
  n_org  <- ncol(mat_org)
  n_tiss <- ncol(mat_tiss)

  cor_mat  <- matrix(NA, nrow = n_org, ncol = n_tiss,
                     dimnames = list(colnames(mat_org), colnames(mat_tiss)))
  pval_mat <- matrix(NA, nrow = n_org, ncol = n_tiss,
                     dimnames = list(colnames(mat_org), colnames(mat_tiss)))
  n_mat    <- get_pairwise_n_matrix(mat_org, mat_tiss)

  for (i in seq_len(n_org)) {
    for (j in seq_len(n_tiss)) {
      x    <- mat_org[, i]
      y    <- mat_tiss[, j]
      keep <- !is.na(x) & !is.na(y)
      n    <- sum(keep)

      if (n < min_n) {
        log_info(sprintf(
          "Masking pair [%s x %s]: only %d shared genes (min_n=%d)",
          colnames(mat_org)[i], colnames(mat_tiss)[j], n, min_n
        ))
        next
      }

      test_result    <- cor.test(x[keep], y[keep], method = method)
      cor_mat[i, j]  <- test_result$estimate
      pval_mat[i, j] <- test_result$p.value
    }
  }

  # FDR correction across all pairs
  pval_adj_mat <- matrix(
    p.adjust(as.vector(pval_mat), method = "BH"),
    nrow = n_org, ncol = n_tiss,
    dimnames = list(colnames(mat_org), colnames(mat_tiss))
  )

  list(
    cor  = cor_mat,
    pval = pval_mat,
    padj = pval_adj_mat,
    n    = n_mat
  )
}

# ─────────────────────────────────────────────
# 5. VISUALIZATION
# ─────────────────────────────────────────────
plot_correlation_heatmap <- function(
    cor_mat,
    pval_mat,
    n_mat,
    title       = "Organoid vs Tissue",
    subtitle    = "",
    sig_only    = FALSE,
    padj_thresh = 0.05
) {
  # Replace NA with 0 only for clustering
  mat_for_clustering <- cor_mat
  mat_for_clustering[is.na(mat_for_clustering)] <- 0

  display_mat <- cor_mat
  if (sig_only) {
    display_mat[pval_mat > padj_thresh] <- NA
  }

  # Significance stars
  sig_stars <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                      dimnames = dimnames(pval_mat))
  sig_stars[!is.na(pval_mat) & pval_mat < 0.05]  <- "*"
  sig_stars[!is.na(pval_mat) & pval_mat < 0.01]  <- "**"
  sig_stars[!is.na(pval_mat) & pval_mat < 0.001] <- "***"

  # Cell annotation: r + stars
  cell_labels <- matrix("NA", nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                        dimnames = dimnames(cor_mat))
  not_na <- !is.na(cor_mat)
  cell_labels[not_na] <- sprintf("%.2f%s", cor_mat[not_na], sig_stars[not_na])

  limit <- max(abs(cor_mat), na.rm = TRUE)

  # Pre-compute clustering order on NA-imputed matrix
  row_order <- hclust(dist(mat_for_clustering))$order
  col_order <- hclust(dist(t(mat_for_clustering)))$order

  full_title <- if (nchar(subtitle) > 0) paste0(title, "\n", subtitle) else title

  # --- Plot 1: Correlation heatmap ---
  pheatmap(
    display_mat[row_order, col_order],
    display_numbers = cell_labels[row_order, col_order],
    number_color    = "black",
    color           = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
    breaks          = seq(-limit, limit, length.out = 101),
    border_color    = "grey90",
    na_col          = "grey85",
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    fontsize        = 11,
    fontsize_number = 8,
    main            = full_title,
    legend_breaks   = c(-limit, -limit/2, 0, limit/2, limit),
    legend_labels   = c(
      sprintf("%.2f", -limit),
      sprintf("%.2f", -limit/2),
      "0",
      sprintf("%.2f",  limit/2),
      sprintf("%.2f",  limit)
    ),
    angle_col = 45
  )

  # --- Plot 2: N genes per pair ---
  n_display <- n_mat[row_order, col_order]
  n_labels  <- matrix(
    as.character(n_display),
    nrow = nrow(n_display), ncol = ncol(n_display),
    dimnames = dimnames(n_display)
  )

  pheatmap(
    n_display,
    display_numbers = n_labels,
    number_color    = "black",
    color           = colorRampPalette(c("#f7fbff", "#2171b5"))(100),
    border_color    = "grey90",
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    fontsize        = 11,
    fontsize_number = 8,
    main            = paste0("N genes per pair\n", subtitle),
    angle_col       = 45
  )
}
# ─────────────────────────────────────────────
# 6. QC REPORT
# ─────────────────────────────────────────────

correlation_qc_report <- function(cor_result) {
  cat("\n", strrep("=", 60), "\n")
  cat("  CORRELATION QC REPORT\n")
  cat(strrep("=", 60), "\n\n")

  cor_vals <- as.vector(cor_result$cor)
  cor_vals <- cor_vals[!is.na(cor_vals)]

  cat(sprintf("  Pairs computed:        %d of %d\n",
              sum(!is.na(cor_result$cor)), length(cor_result$cor)))
  cat(sprintf("  Pairs masked (low N):  %d\n",
              sum(is.na(cor_result$cor))))
  cat(sprintf("  Mean r:                %.3f\n", mean(cor_vals)))
  cat(sprintf("  Std r:                 %.3f\n", sd(cor_vals)))
  cat(sprintf("  Min r:                 %.3f\n", min(cor_vals)))
  cat(sprintf("  Max r:                 %.3f\n", max(cor_vals)))
  cat(sprintf("  Sig pairs (padj<0.05): %d of %d (%.1f%%)\n",
              sum(cor_result$padj < 0.05, na.rm = TRUE),
              sum(!is.na(cor_result$padj)),
              100 * mean(cor_result$padj < 0.05, na.rm = TRUE)))

  cat("\n[Best Match per Organoid Cluster — Spearman]\n")
  best_df <- data.frame(
    organoid    = rownames(cor_result$cor),
    best_tissue = colnames(cor_result$cor)[
      apply(cor_result$cor, 1, which.max)
    ],
    max_r = round(apply(cor_result$cor, 1, max, na.rm = TRUE), 3),
    min_r = round(apply(cor_result$cor, 1, min, na.rm = TRUE), 3)
  ) %>%
    mutate(delta = round(max_r - min_r, 3)) %>%
    arrange(desc(max_r))

  print(best_df, row.names = FALSE)
  cat("\n  [delta = discriminability: higher = cleaner mapping]\n")
  cat(strrep("=", 60), "\n\n")
}

# ─────────────────────────────────────────────
# 7. RUN
# ─────────────────────────────────────────────

CHOSEN_METHOD    <- "union"
CHOSEN_STAT      <- "t"
USE_CONC         <- FALSE
MIN_N_PAIRS      <- 10
USE_ADJ_ORG      <- FALSE   # nominal — organoid underpowered
USE_ADJ_TISS     <- FALSE   # nominal — tissue BD effect too small for FDR
PVAL_THRES_ORG   <- 0.05
PVAL_THRES_TISS  <- 0.05
USE_SIG_FILTER   <- TRUE    # TRUE = sig filter, FALSE = all shared genes

org_df_flat  <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(
  org_df_flat, tiss_df_flat,
  method          = CHOSEN_METHOD,
  stat_type       = CHOSEN_STAT,
  use_concordance = USE_CONC,
  pval_thres_org  = PVAL_THRES_ORG,
  pval_thres_tiss = PVAL_THRES_TISS,
  use_adj_org     = USE_ADJ_ORG,
  use_adj_tiss    = USE_ADJ_TISS,
  use_sig_filter  = USE_SIG_FILTER,
  min_n_pairs     = MIN_N_PAIRS
)
log_info(sprintf("Final features: %d genes", length(final_features)))

mat_org  <- prep_cluster_matrix(org_df_flat,  final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes   <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final  <- as.matrix(mat_org[common_genes,  , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

log_info(sprintf("Matrix dims — Organoid: %s | Tissue: %s",
                 paste(dim(mat_org_final),  collapse = "x"),
                 paste(dim(mat_tiss_final), collapse = "x")))

cor_result <- compute_correlation_with_pval(
  mat_org_final, mat_tiss_final,
  method = "spearman",
  min_n  = MIN_N_PAIRS
)

correlation_qc_report(cor_result)

plot_subtitle <- sprintf(
  "Method: %s | Stat: %s | Concordance: %s | Genes: %d\nOrganoid: %s p<%.2f | Tissue: %s p<%.2f",
  toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),
  ifelse(USE_ADJ_ORG,  "FDR", "nominal"), PVAL_THRES_ORG,
  ifelse(USE_ADJ_TISS, "FDR", "nominal"), PVAL_THRES_TISS
)

plot_correlation_heatmap(
  cor_mat     = cor_result$cor,
  pval_mat    = cor_result$padj,
  n_mat       = cor_result$n,
  title       = "Organoid vs Tissue — Spearman Correlation",
  subtitle    = plot_subtitle,
  sig_only    = FALSE,
  padj_thresh = 0.05
)

# Confirm sign flip by checking a known BD marker
# BDNF, SCN1A, CACNA1C are commonly reported in BD
known_bd_genes <- c("BDNF", "CACNA1C", "SCN1A", "ANK3", "NRXN1")

cat("=== Known BD genes direction check ===\n")
for (gene in known_bd_genes) {
  if (gene %in% rownames(mat_org_final) & gene %in% rownames(mat_tiss_final)) {
    cat(sprintf("\n%s:\n", gene))
    cat("  Organoid t:"); print(round(mat_org_final[gene, ], 3))
    cat("  Tissue t:  "); print(round(mat_tiss_final[gene, ], 3))
  }
}

# Check original dreamlet result column names
cat("\n=== Organoid df columns ===\n")
print(colnames(org_df_flat))

cat("\n=== Tissue df columns ===\n")
print(colnames(tiss_df_flat))

# Check what the coef column says
cat("\n=== Organoid unique coef values ===\n")
if ("coef" %in% colnames(org_df_flat))  print(unique(org_df_flat$coef))
if ("L1"   %in% colnames(org_df_flat))  print(unique(org_df_flat$L1))
if ("term"  %in% colnames(org_df_flat)) print(unique(org_df_flat$term))

cat("\n=== Tissue unique coef values ===\n")
if ("coef" %in% colnames(tiss_df_flat))  print(unique(tiss_df_flat$coef))
if ("L1"   %in% colnames(tiss_df_flat))  print(unique(tiss_df_flat$L1))
if ("term"  %in% colnames(tiss_df_flat)) print(unique(tiss_df_flat$term))

# Quick comparison: t vs z.std correlation for Astro vs Astro
for (stat in c("t", "z.std", "logFC")) {
  mat_o <- prep_cluster_matrix(org_df_flat,  final_features, stat)
  mat_ti <- prep_cluster_matrix(tiss_df_flat, final_features, stat)
  cg <- sort(intersect(rownames(mat_o), rownames(mat_ti)))
  x  <- as.matrix(mat_o[cg,  , drop=FALSE])[, "Astro"]
  y  <- as.matrix(mat_ti[cg, , drop=FALSE])[, "Astro"]
  keep <- !is.na(x) & !is.na(y)
  cat(sprintf("%-8s — Astro vs Astro r=%.4f (n=%d)\n",
              stat,
              cor(x[keep], y[keep], method="spearman"),
              sum(keep)))
}

# Is r=0.116 with n=4630 actually significant?
cor.test(
  as.numeric(mat_org_final[, "Astro"]),
  as.numeric(mat_tiss_final[, "Astro"]),
  method = "spearman"
)


library(fgsea)

# ─────────────────────────────────────────────
# Load local GMT files
# ─────────────────────────────────────────────

load_gmt <- function(gmt_path) {
  lines    <- readLines(gmt_path)
  pathways <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    genes <- parts[3:length(parts)]
    genes[genes != ""]
  })
  names(pathways) <- sapply(lines, function(line) {
    strsplit(line, "\t")[[1]][1]
  })
  return(pathways)
}

gmt_gobp <- load_gmt("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c5.go.bp.v7.5.1.symbols.gmt")
gmt_kegg <- load_gmt("../GMT_Files/msigdb_v7.5.1_files_to_download_locally/msigdb_v7.5.1_GMTs/c2.cp.kegg.v7.5.1.symbols.gmt")

# ─────────────────────────────────────────────
# ORA — hypergeometric test
# ─────────────────────────────────────────────

run_ora <- function(
    gene_list,
    pathways,
    universe,
    min_size    = 5,
    max_size    = 500,
    padj_thresh = 0.05
) {
  n_universe <- length(universe)
  n_query    <- length(gene_list)

  results <- map_dfr(names(pathways), function(pw_name) {
    pw_genes  <- intersect(pathways[[pw_name]], universe)
    n_pw      <- length(pw_genes)

    if (n_pw < min_size | n_pw > max_size) return(NULL)

    overlap   <- intersect(gene_list, pw_genes)
    n_overlap <- length(overlap)

    if (n_overlap == 0) return(NULL)

    pval <- phyper(
      n_overlap - 1,
      n_pw,
      n_universe - n_pw,
      n_query,
      lower.tail = FALSE
    )

    data.frame(
      Term      = pw_name,
      n_pathway = n_pw,
      n_overlap = n_overlap,
      n_query   = n_query,
      pval      = pval,
      Genes     = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    filter(padj < padj_thresh) %>%
    arrange(padj)

  return(results)
}

# ─────────────────────────────────────────────
# Print helper — explicit dplyr:: to avoid
# AnnotationDbi::select namespace collision
# ─────────────────────────────────────────────

print_ora <- function(ora_result, n = 10) {
  if (is.null(ora_result) || nrow(ora_result) == 0) {
    cat("  No significant terms found.\n")
    return(invisible(NULL))
  }
  ora_result %>%
    dplyr::select(Term, n_pathway, n_overlap, padj, Genes) %>%
    head(n) %>%
    print(row.names = FALSE)
}

# ─────────────────────────────────────────────
# Run ORA
# ─────────────────────────────────────────────

universe <- intersect(unique(org_df_flat$ID), unique(tiss_df_flat$ID))
cat(sprintf("Universe size: %d genes\n", length(universe)))

cat("\n=== ORA: Astro concordant — GO BP ===\n")
ora_astro_gobp <- run_ora(
  gene_list = concordant_results$Astro_Astro$gene,
  pathways  = gmt_gobp,
  universe  = universe
)
print_ora(ora_astro_gobp)

cat("\n=== ORA: Astro concordant — KEGG ===\n")
ora_astro_kegg <- run_ora(
  gene_list = concordant_results$Astro_Astro$gene,
  pathways  = gmt_kegg,
  universe  = universe
)
print_ora(ora_astro_kegg)

cat("\n=== ORA: Neuron/IN concordant — GO BP ===\n")
ora_neuron_gobp <- run_ora(
  gene_list = concordant_results$Neuron_IN$gene,
  pathways  = gmt_gobp,
  universe  = universe
)
print_ora(ora_neuron_gobp)

cat("\n=== ORA: Neuron/IN concordant — KEGG ===\n")
ora_neuron_kegg <- run_ora(
  gene_list = concordant_results$Neuron_IN$gene,
  pathways  = gmt_kegg,
  universe  = universe
)
print_ora(ora_neuron_kegg)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

# --- 1. UTILITY & PROCESSING FUNCTIONS ---

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
  flush(stderr())
}

#' Helper to extract gene lists from Bioconductor DataFrames
#' @param df S4 DataFrame or Tibble
#' @param prefix String to prepend to cell type names
get_gene_list <- function(df, prefix) {
    if(prefix=='O_'){
      df %>%
        as.data.frame() %>% # Cast S4 DataFrame to base R data.frame
        as_tibble() %>%
        filter(P.Value < 0.05& abs(logFC) > 0) %>%
        mutate(set_name = paste0(prefix, assay)) %>%
        split(.$set_name) %>%
        map(~ .x$ID)
    }else{
       df %>%
            as.data.frame() %>% # Cast S4 DataFrame to base R data.frame
            as_tibble() %>%
            filter(adj.P.Val < 0.05& abs(logFC) > 0) %>%
            mutate(set_name = paste0(prefix,assay)) %>%
            split(.$set_name) %>%
            map(~ .x$ID)
      }
          
}

#' Feature Selection: Background, Significance, and Concordance
get_feature_genes <- function(df_org, df_tiss, method = "union", stat_type = "t", use_concordance = TRUE,pval_thres=0.05) {
  # 1. Background: Genes physically present in both
  common_background <- intersect(unique(df_org$ID), unique(df_tiss$ID))
  
  # 2. Significance: Nominal P < 0.05
  sig_ids_org <- unique(df_org$ID[df_org$P.Value < pval_thres])
  sig_ids_tiss <- unique(df_tiss$ID[df_tiss$P.Value < pval_thres])
  
  if (method == "union") {
    target_sig <- union(sig_ids_org, sig_ids_tiss)
  } else {
    target_sig <- intersect(sig_ids_org, sig_ids_tiss)
  }
  
  candidate_genes <- intersect(target_sig, common_background)
  
  # 3. Optional Concordance: Same direction based on chosen stat_type
  if (use_concordance) {
    trend_org <- df_org %>% filter(ID %in% candidate_genes) %>% 
      group_by(ID) %>% summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop")
    
    trend_tiss <- df_tiss %>% filter(ID %in% candidate_genes) %>% 
      group_by(ID) %>% summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop")
    
    candidate_genes <- inner_join(trend_org, trend_tiss, by = "ID") %>%
      filter(sign(val.x) == sign(val.y)) %>%
      pull(ID)
  }
  
  return(candidate_genes)
}

#' Matrix Preparation: Pivot to Genes x Clusters using chosen stat_type
prep_cluster_matrix <- function(df, genes, stat_type) {
  df %>%
    as.data.frame() %>%
    filter(ID %in% genes) %>%
    group_by(assay, ID) %>%
    summarize(val = mean(get(stat_type), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = assay, values_from = val) %>%
    column_to_rownames("ID")
}

get_pairwise_n_matrix <- function(mat_x, mat_y) {
  indicator_x <- !is.na(mat_x); indicator_y <- !is.na(mat_y)
  return(as.matrix(t(indicator_x) %*% indicator_y))
}

COEF <- "Diff"
organoid_df <- topTable(readRDS('2026_03_06_Scanpy_020_2dpass_annotated_CLASS_cova_DiffPDvsCTRL_1.9.rds'), coef = COEF, number = Inf)
tissue_df <- topTable(readRDS('Lyra_data/dreamlet_class_upitt.rds'), coef = COEF, number = Inf)
#tissue_df <- readRDS('Lyra_data/subclass_upitt_dreamlet_toptable.rds')
tissue_sets <- get_gene_list(tissue_df, "T_")
organoid_sets <- get_gene_list(organoid_df, "O_")


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "t"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.05         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s \n Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)

options(repr.plot.width = 8, repr.plot.height = 5, repr.plot.res = 200)

pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize_number = 8,
  fontsize_main = 10
)
pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "logFC"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.05         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s \n Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)

options(repr.plot.width = 8, repr.plot.height = 5, repr.plot.res = 200)

pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize_number = 8,
  fontsize_main = 10
)
pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)

library(dplyr)
library(fgsea)
library(ggplot2)
library(stringr)

# ---------------------------------------------------------
# STEP 3: Format Ranked Lists for ALL Clusters
# ---------------------------------------------------------
cat(sprintf("[%s] INFO: Formatting full ranked gene lists across all clusters...\n", Sys.time()), file=stderr())

# Split the full unthresholded dataframe by cluster, then extract named, sorted vectors
ranks_list <- lapply(split(fast_markers, fast_markers$cluster), function(df) {
  ranks <- df$avg_log2FC
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
})

# CHECK: Ensure the list has slots for all your clusters and contains thousands of genes
length(ranks_list)
names(ranks_list)
length(ranks_list[[1]]) # Should be ~10,000+ genes, not just the top hits


# ---------------------------------------------------------
# STEP 4: Run Multi-Cluster fgsea
# ---------------------------------------------------------
cat(sprintf("[%s] INFO: Running fgseaMultilevel across %d clusters...\n", Sys.time(), length(ranks_list)), file=stderr())
# ---------------------------------------------------------
# THE FIXED FGSEA CALL
# ---------------------------------------------------------
cat(sprintf("[%s] INFO: Running fgseaMultilevel with positive-only optimized math...\n", Sys.time()), file=stderr())

fgsea_all_list <- lapply(names(ranks_list), function(clust) {
  res <- fgseaMultilevel(
    pathways = pathways, 
    stats = ranks_list[[clust]],
    minSize = 15,
    maxSize = 500,
    scoreType = "pos",       # FIX 1: Tells the math to only expect positive logFC values
    nPermSimple = 10000      # FIX 2: Increases permutations to stabilize extreme p-values
  )
  
  res$Cluster <- clust
  return(res)
})


# Bind everything into one master dataframe and filter for significance
fgsea_all_df <- dplyr::bind_rows(fgsea_all_list) %>%
  dplyr::filter(padj < 0.05)

# CHECK: View the top hits across the object
head(fgsea_all_df[, c("Cluster", "pathway", "NES", "padj")])


# ---------------------------------------------------------
# STEP 5: Plot the Enrichment per Cluster (Custom DotPlot)
# ---------------------------------------------------------
cat(sprintf("[%s] INFO: Generating multi-cluster fgsea DotPlot...\n", Sys.time()), file=stderr())

# Extract the top 5 enriched pathways per cluster by NES (Normalized Enrichment Score)
top_fgsea <- fgsea_all_df %>%
  dplyr::group_by(Cluster) %>%
  dplyr::slice_max(order_by = NES, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    # Clean up pathway names for the plot
    clean_pathway = stringr::str_replace_all(pathway, "^[A-Z0-9]+_", ""),
    clean_pathway = stringr::str_trunc(clean_pathway, width = 50)
  )



# # Build the DotPlot natively in ggplot2
# final_plot <- ggplot(top_fgsea, aes(x = Cluster, y = clean_pathway)) +
#   # Size maps to Enrichment Score, Color maps to Significance
#   geom_point(aes(size = NES, color = -log10(padj))) +
#   scale_color_viridis_c(option = "plasma") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
#   ) +
#   labs(
#     title = "Top 5 Enriched Cell Types per Cluster (fgsea)",
#     x = "Organoid Cluster",
#     y = "MSigDB C8 Cell Type Signature",
#     size = "NES",
#     color = "-log10(padj)"
#   )

# # View and Save
# # print(final_plot)
# # ggsave("all_clusters_fgsea_dotplot.pdf", plot = final_plot, width = 14, height = 10)

# options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 300)
# print(final_plot)

library(tidyr)
library(tibble)

cat(sprintf("[%s] INFO: Performing hierarchical clustering for plot ordering...\n", Sys.time()), file=stderr())

# 1. Pivot to a wide matrix of NES scores (Pathways = Rows, Clusters = Columns)
# We fill missing enrichment values with 0 so the math works
nes_matrix <- top_fgsea %>%
  dplyr::select(Cluster, clean_pathway, NES) %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = NES, values_fill = list(NES = 0)) %>%
  tibble::column_to_rownames("clean_pathway") %>%
  as.matrix()

# 2. Mathematically cluster the Pathways (Rows)
pathway_dist <- dist(nes_matrix, method = "euclidean")
pathway_hclust <- hclust(pathway_dist, method = "ward.D2")
pathway_order <- rownames(nes_matrix)[pathway_hclust$order]

# 3. Mathematically cluster the Organoid Clusters (Columns)
cluster_dist <- dist(t(nes_matrix), method = "euclidean")
cluster_hclust <- hclust(cluster_dist, method = "ward.D2")
cluster_order <- colnames(nes_matrix)[cluster_hclust$order]

# 4. Lock in the new clustered order by updating the factor levels
top_fgsea_clustered <- top_fgsea %>%
  dplyr::mutate(
    clean_pathway = factor(clean_pathway, levels = pathway_order),
    Cluster = factor(Cluster, levels = cluster_order)
  )

cat(sprintf("[%s] INFO: Generating hierarchically clustered DotPlot...\n", Sys.time()), file=stderr())



# 5. Plot using the new ordered dataframe
final_plot_clustered <- ggplot(top_fgsea_clustered, aes(x = Cluster, y = clean_pathway)) +
  geom_point(aes(size = NES, color = -log10(padj))) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    title = "Top 5 Enriched Cell Types per Cluster (Hierarchically Clustered)",
    x = "Organoid Cluster (Grouped by Similarity)",
    y = "MSigDB C8 Cell Type Signature",
    size = "NES",
    color = "-log10(padj)"
  )

# print(final_plot_clustered)

print(final_plot_clustered)

options(repr.plot.width = 10, repr.plot.height = 15, repr.plot.res = 300)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

# 1. Filter for PanglaoDB results
panglao_results <- listofresults[grep("PanglaoDB", names(listofresults))]

# 2. Extract and Process Data
plot_data <- imap_dfr(panglao_results, function(res, name) {
  cluster_id <- str_extract(name, "\\d+$")
  
  # Convert to dataframe
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  
  # Robust column detection: find the column that looks like adjusted p-value
  # This handles 'p.adjust', 'Adjusted.P.value', 'padj', etc.
  p_col <- grep("adj|p.*val", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  term_col <- grep("Term|Description", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  ratio_col <- grep("Ratio|Overlap", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  # If we can't find a p-value column, skip this cluster
  if (is.na(p_col)) return(NULL)
  
  df %>%
    # Rename columns to standard names for plotting
    rename(p_val = !!p_col, Term = !!term_col) %>%
    # Use GeneRatio if it exists, otherwise use a placeholder for size
    mutate(Ratio = if(!is.na(ratio_col)) as.numeric(str_extract(get(ratio_col), "^\\d+")) else 1) %>%
    slice_min(order_by = p_val, n = 5, with_ties = FALSE) %>%
    mutate(Cluster = as.numeric(cluster_id))
})

# 3. Create the Combined Dot Plot
dotplot <- ggplot(plot_data, aes(x = factor(Cluster), y = reorder(Term, -p_val))) +
  geom_point(aes(size = Ratio, color = p_val)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "Top 5 PanglaoDB Terms per Cluster",
    x = "Cluster ID",
    y = "Enriched Cell Types",
    size = "Enrichment Score/Ratio",
    color = "Adj. P-value"
  ) +
  theme(axis.text.y = element_text(size = 8))

# Logging for reproducibility
cat(sprintf("[%s] Pipeline: Processed %d clusters.\n", Sys.time(), length(unique(plot_data$Cluster))), file=stderr())

print(dotplot)

listofresults <- readRDS('../listofresults.rds')
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 300)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)

# 1. Filter for PanglaoDB results
panglao_results <- listofresults[grep("PanglaoDB", names(listofresults))]

# 1. Process Data
plot_data <- imap_dfr(panglao_results, function(res, name) {
  cluster_id <- str_extract(name, "\\d+$")
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  
  p_col <- grep("adj|p.*val", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  term_col <- grep("Term|Description", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  ratio_col <- grep("Ratio|Overlap", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  df %>%
    rename(p_val = !!p_col, Term = !!term_col) %>%
    mutate(Ratio = if(!is.na(ratio_col)) as.numeric(str_extract(get(ratio_col), "^\\d+")) else 1) %>%
    slice_min(order_by = p_val, n = 5, with_ties = FALSE) %>%
    mutate(Cluster = as.numeric(cluster_id))
})

# 2. Hierarchical Clustering of Clusters (X-axis)
# Using explicit namespaces to avoid AnnotationDbi conflicts
cluster_matrix <- plot_data %>%
  dplyr::select(Term, Cluster) %>% 
  mutate(Present = 1) %>%
  pivot_wider(names_from = Cluster, values_from = Present, values_fill = 0) %>%
  tibble::column_to_rownames("Term")

# Compute distance and reorder
cluster_dist <- dist(t(cluster_matrix)) 
cluster_hclust <- hclust(cluster_dist)
cluster_order <- cluster_hclust$labels[cluster_hclust$order]

# Apply the order
plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)

# 3. Plotting
dotplot <- ggplot(plot_data, aes(x = Cluster, y = reorder(Term, -p_val))) +
  geom_point(aes(size = Ratio, color = p_val)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "PanglaoDB: Clusters Grouped by Similarity",
    x = "Cluster ID (Clustered)",
    y = "Enriched Terms",
    size = "Enrichment Score",
    color = "Adj. P-value"
  ) +
  theme(axis.text.x = element_text(face = "bold"))

cat(sprintf("[%s] SUCCESS | Clusters reordered by hierarchical similarity.\n", Sys.time()), file=stderr())

print(dotplot)

names(listofresults)


fast_markers10<-fast_markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) #100  gives best results

fast_markers10

listofresults <- readRDS('../Brain_listofresults.rds')
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)

# 1. Filter for PanglaoDB results
panglao_results <- listofresults[grep("Azimuth_2023", names(listofresults))]

# 1. Process Data
# 1. Process Data (Updated with Filter)
plot_data <- imap_dfr(panglao_results, function(res, name) {
  cluster_id <- str_extract(name, "\\d+$")
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  
  p_col <- grep("adj|p.*val", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  term_col <- grep("Term|Description", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  ratio_col <- grep("Ratio|Overlap", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  df %>%
    rename(p_val = !!p_col, Term = !!term_col) %>%
    # Remove terms containing "Mouse" (case-insensitive)
    filter(!str_detect(Term, regex("Mouse|down", ignore_case = TRUE))) %>%
    mutate(Ratio = if(!is.na(ratio_col)) as.numeric(str_extract(get(ratio_col), "^\\d+")) else 1) %>%
    # Recalculate top 5 AFTER filtering
    slice_min(order_by = p_val, n = 5, with_ties = FALSE) %>%
    mutate(Cluster = as.numeric(cluster_id))
})

# 2. Hierarchical Clustering (Refined Matrix)
cluster_matrix <- plot_data %>%
  dplyr::select(Term, Cluster) %>% 
  mutate(Present = 1) %>%
  pivot_wider(
    names_from = Cluster, 
    values_from = Present, 
    values_fill = 0,
    values_fn = list(Present = length) 
  ) %>%
  mutate(across(-Term, ~ as.numeric(.x > 0))) %>%
  tibble::column_to_rownames("Term")

cluster_dist <- dist(t(cluster_matrix)) 
cluster_hclust <- hclust(cluster_dist)
cluster_order <- cluster_hclust$labels[cluster_hclust$order]

plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)

# 3. Plotting
dotplot <- ggplot(plot_data, aes(x = Cluster, y = reorder(Term, -p_val))) +
  geom_point(aes(size = Ratio, color = p_val)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "Allen Brain Atlas: Top 5 (Non-Mouse) Terms",
    subtitle = "Hierarchical clustering of clusters based on shared regional/cell markers",
    x = "Cluster ID (Clustered)",
    y = "Enriched Brain Regions/Cells",
    size = "Overlap/Ratio",
    color = "Adj. P-value"
  ) +
  theme(axis.text.x = element_text(face = "bold"))

print(dotplot)

listofresults <- readRDS('../Brain_listofresults.rds')
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)

# 1. Filter for PanglaoDB results
panglao_results <- listofresults[grep("PanglaoDB_Augmented_2021", names(listofresults))]

# 1. Process Data
# 1. Process Data (Updated with Filter)
plot_data <- imap_dfr(panglao_results, function(res, name) {
  cluster_id <- str_extract(name, "\\d+$")
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  
  p_col <- grep("adj|p.*val", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  term_col <- grep("Term|Description", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  ratio_col <- grep("Ratio|Overlap", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  df %>%
    rename(p_val = !!p_col, Term = !!term_col) %>%
    # Remove terms containing "Mouse" (case-insensitive)
    filter(!str_detect(Term, regex("Mouse|down", ignore_case = TRUE))) %>%
    mutate(Ratio = if(!is.na(ratio_col)) as.numeric(str_extract(get(ratio_col), "^\\d+")) else 1) %>%
    # Recalculate top 5 AFTER filtering
    slice_min(order_by = p_val, n = 5, with_ties = FALSE) %>%
    mutate(Cluster = as.numeric(cluster_id))
})

# 2. Hierarchical Clustering (Refined Matrix)
cluster_matrix <- plot_data %>%
  dplyr::select(Term, Cluster) %>% 
  mutate(Present = 1) %>%
  pivot_wider(
    names_from = Cluster, 
    values_from = Present, 
    values_fill = 0,
    values_fn = list(Present = length) 
  ) %>%
  mutate(across(-Term, ~ as.numeric(.x > 0))) %>%
  tibble::column_to_rownames("Term")

cluster_dist <- dist(t(cluster_matrix)) 
cluster_hclust <- hclust(cluster_dist)
cluster_order <- cluster_hclust$labels[cluster_hclust$order]

plot_data$Cluster <- factor(plot_data$Cluster, levels = cluster_order)

# 3. Plotting
dotplot <- ggplot(plot_data, aes(x = Cluster, y = reorder(Term, -p_val))) +
  geom_point(aes(size = Ratio, color = p_val)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "Allen Brain Atlas: Top 5 (Non-Mouse) Terms",
    subtitle = "Hierarchical clustering of clusters based on shared regional/cell markers",
    x = "Cluster ID (Clustered)",
    y = "Enriched Brain Regions/Cells",
    size = "Overlap/Ratio",
    color = "Adj. P-value"
  ) +
  theme(axis.text.x = element_text(face = "bold"))

print(dotplot)

library(SingleR)
pred.hesc <- SingleR(test = sce, labels = sce$leiden_08_set)

sce <- as.SingleCellExperiment(seurat_obj)
pb <- aggregateToPseudoBulk(sce,
  assay = "counts",
  cluster_id = "leiden_08_donor",
  sample_id = "final_donor",
  verbose = FALSE
)
df_cts <- cellTypeSpecificity(pb)

df_cts <- cellTypeSpecificity(pb)
# retain only genes with total CPM summed across cell type > 100
df_cts <- df_cts[df_cts$totalCPM > 100, ]
# Define N, the number of top genes you want (in this case, 5)
N <- 500

# Use apply() with MARGIN=2 (columns) and a custom function
# The function sorts the column and returns the top N row names (genes)
top_n_genes <- apply(df_cts, 2, function(col_values) {
    # 1. Use 'order' with 'decreasing=TRUE' to get the indices 
    #    of the values sorted from largest to smallest.
    # 2. Select the first N indices (the top N values).
    top_indices <- order(col_values, decreasing = TRUE)[1:N]
    
    # 3. Use these indices to select the corresponding row names (genes)
    #    from the full data frame/matrix's row names.
    return(rownames(df_cts)[top_indices])
})

df_long <- as.data.frame(top_n_genes) %>%
  dplyr::select(-1) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Cluster",
    values_to = "genes"
  ) %>%
  filter(!is.na(genes) & genes != "") %>%
  # Ensure Cluster is treated as a factor or numeric for correct sorting
  mutate(Cluster = factor(Cluster, levels = unique(Cluster)))

for_enrichR <- as.data.frame(df_long)
colnames(for_enrichR) <- c('cluster','gene')

head(for_enrichR)
saveRDS(for_enrichR,'../for_enrichR.rds')

#rm(seurat_obj)

top_n_genes

# res_dl_dx1 <- readRDS('CLASS_cova_DiffPDvsCTRL_1.9.rds')
res_dl_dx1 <- readRDS('2026_03_03_Scanpy_020_2dpass_annotated_CLASS_cova_DiffPDvsCTRL_1.9.rds')

summary_tab <- topTable(res_dl_dx1, coef = COEF, number = Inf) %>%
  as_tibble() %>%
  group_by(assay) %>%
  summarize(
    nDE = sum(P.Value < 0.05),
    pi1 = 1 - pi0est(P.Value)$pi0,
    nGenes = length(adj.P.Val)
  ) %>%
  mutate(assay = factor(assay, assayNames(res_dl_dx1)))
summary_tab

library(zellkonverter)
library(SingleCellExperiment)
library(dreamlet)
library(ggplot2)
library(tidyverse)
library(aplot)
library(ggtree)
library(scattermore)
library(zenith)
library(crumblr)
library(GSEABase)
library(qvalue)
library(BiocParallel)
library(cowplot)
library(DelayedArray)
library(tidyverse)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
# --- Print Versions for Reproducibility ---

cat("Session package versions:\n")
cat(paste0("dreamlet v", packageVersion("dreamlet")), "\n")
cat(paste0("crumblr v", packageVersion("crumblr")), "\n")
cat(paste0("variancePartition v", packageVersion("variancePartition")), "\n")
cat(paste0("zenith v", packageVersion("zenith")), "\n")
cat(paste0("zellkonverter v", packageVersion("zellkonverter")), "\n")
cat(paste0("BiocManager v", BiocManager::version()), "\n")

print(getwd())
library(Seurat)
options(Seurat.object.assay.version = 'v3')
setwd(OUTDIR)
# --- Step 1: Load H5AD and Metadata ---

#' Helper to extract gene lists from Bioconductor DataFrames
#' @param df S4 DataFrame or Tibble
#' @param prefix String to prepend to cell type names
get_gene_list <- function(df, prefix) {
    if(prefix=='O_'){
      df %>%
        as.data.frame() %>% # Cast S4 DataFrame to base R data.frame
        as_tibble() %>%
        filter(P.Value < 0.05& abs(logFC) > 0) %>%
        mutate(set_name = paste0(prefix, assay)) %>%
        split(.$set_name) %>%
        map(~ .x$ID)
    }else{
       df %>%
            as.data.frame() %>% # Cast S4 DataFrame to base R data.frame
            as_tibble() %>%
            filter(adj.P.Val < 0.05& abs(logFC) > 0) %>%
            mutate(set_name = paste0(prefix,assay)) %>%
            split(.$set_name) %>%
            map(~ .x$ID)
      }
          
}

# --- FIX: Plotting function for directional counts with Integer Legends ---
plt_overalap_degs <- function(res_mats, nominal_a = TRUE, nominal_b = TRUE,
                              plot.width = 12, plot.height = 4, plot.res = 200) {
    max_up <- max(res_mats$up_counts)
    max_dn <- max(res_mats$dn_counts)
    
    col_up <- colorRamp2(c(0, max_up), c("white", "#D73027"))
    col_dn <- colorRamp2(c(0, max_dn), c("white", "#4575B4"))
    
    # Legend with rounded/integer labels
    at_up <- unique(round(seq(0, max_up, length.out = 3)))
    at_dn <- unique(round(seq(0, max_dn, length.out = 3)))
    
    lgd_up = Legend(title = "# Up-reg", col_fun = col_up, at = at_up, labels = as.integer(at_up))
    lgd_dn = Legend(title = "# Down-reg", col_fun = col_dn, at = at_dn, labels = as.integer(at_dn))
    
    column_title <- if (nominal_a) "Tissue subclass (p.val)" else "Tissue subclass (p.adj)"
    row_title <- if (nominal_b) "Organoid subclass (p.val)" else "Organoid subclass (p.adj)"
        
    ht <- Heatmap(res_mats$up_counts, 
                  border = TRUE,
                  rect_gp = gpar(type = "none"),
                  column_title = column_title, row_title = row_title,
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  show_heatmap_legend = FALSE,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      c_up <- res_mats$up_counts[i, j]
                      c_dn <- res_mats$dn_counts[i, j]
                      grid.polygon(x = c(x - width/2, x + width/2, x + width/2),
                                   y = c(y + height/2, y + height/2, y - height/2),
                                   gp = gpar(fill = col_up(c_up), col = "white", lwd = 0.5))
                      grid.polygon(x = c(x - width/2, x + width/2, x - width/2),
                                   y = c(y + height/2, y - height/2, y - height/2),
                                   gp = gpar(fill = col_dn(c_dn), col = "white", lwd = 0.5))
                      if(c_up > 0) grid.text(c_up, x + width/4, y + height/4, gp = gpar(fontsize = 8, col = ifelse(c_up > (max_up*0.6), "white", "black")))
                      if(c_dn > 0) grid.text(c_dn, x - width/4, y - height/4, gp = gpar(fontsize = 8, col = ifelse(c_dn > (max_dn*0.6), "white", "black")))
                  })
    options(repr.plot.width = plot.width, repr.plot.height = plot.height, repr.plot.res = plot.res)
    draw(ht, annotation_legend_list = list(lgd_up, lgd_dn))
}

plt_total_vs_directional <- function(res_mats, df_a, df_b, nominal_a=TRUE, nominal_b=TRUE,
                              plot.width = 12, plot.height = 4, plot.res = 200) {
    # 1. Calculate Total Overlap (Non-directional)
    p_col_a <- if(nominal_a) "P.Value" else "adj.P.Val"
    p_col_b <- if(nominal_b) "P.Value" else "adj.P.Val"
    assays_a <- rownames(res_mats$up_counts); assays_b <- colnames(res_mats$up_counts)
    
    total_inter_mat <- matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b))
    for(a in assays_a) {
        genes_a <- df_a$ID[df_a$assay == a & df_a[[p_col_a]] < 0.05]
        for(b in assays_b) {
            genes_b <- df_b$ID[df_b$assay == b & df_b[[p_col_b]] < 0.05]
            total_inter_mat[a, b] <- length(intersect(genes_a, genes_b))
        }
    }
    
    # 2. Calculate Directional Agreement (Up-Up + Down-Down)
    agree_mat <- res_mats$up_counts + res_mats$dn_counts
    
    # 3. Setup Plotting
    max_val <- max(total_inter_mat)
    col_total <- colorRamp2(c(0, max_val), c("white", "#6A3D9A")) # Purple for Total
    col_agree <- colorRamp2(c(0, max_val), c("white", "#33A02C")) # Green for Agreement
    
    at_vals <- unique(round(seq(0, max_val, length.out = 3)))
    lgd_tot = Legend(title = "Total Overlap", col_fun = col_total, at = at_vals, labels = as.integer(at_vals))
    lgd_agr = Legend(title = "Concordant DEGs", col_fun = col_agree, at = at_vals, labels = as.integer(at_vals))

    column_title <- if (nominal_a) "Tissue subclass (p.val)" else "Tissue subclass (p.adj)"
    row_title <- if (nominal_b) "Organoid subclass (p.val)" else "Organoid subclass (p.adj)"
        
    ht <- Heatmap(total_inter_mat, 
                  border = TRUE,
                  rect_gp = gpar(type = "none"),
                  column_title = column_title, row_title = row_title,
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  show_heatmap_legend = FALSE,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      v_tot <- total_inter_mat[i, j]
                      v_agr <- agree_mat[i, j]
                      
                      # Top-Right: Total
                      grid.polygon(x = c(x - width/2, x + width/2, x + width/2), 
                                   y = c(y + height/2, y + height/2, y - height/2), 
                                   gp = gpar(fill = col_total(v_tot), col = "white", lwd = 0.5))
                      # Bottom-Left: Agreement
                      grid.polygon(x = c(x - width/2, x + width/2, x - width/2), 
                                   y = c(y + height/2, y - height/2, y - height/2), 
                                   gp = gpar(fill = col_agree(v_agr), col = "white", lwd = 0.5))
                      
                      if(v_tot > 0) grid.text(v_tot, x + width/4, y + height/4, gp = gpar(fontsize = 8, col = ifelse(v_tot > max_val*0.6, "white", "black")))
                      if(v_agr > 0) grid.text(v_agr, x - width/4, y - height/4, gp = gpar(fontsize = 8, col = ifelse(v_agr > max_val*0.6, "white", "black")))
                  })
    
    options(repr.plot.width = plot.width, repr.plot.height = plot.height, repr.plot.res = plot.res)
    draw(ht, annotation_legend_list = list(lgd_tot, lgd_agr))
}


# --- NEW: Non-Directional Total Overlap Counts ---
plt_total_overlap <- function(df_a, df_b, nominal_a=TRUE, nominal_b=TRUE,
                              plot.width = 12, plot.height = 4, plot.res = 200) {
    p_col_a <- if(nominal_a) "P.Value" else "adj.P.Val"
    p_col_b <- if(nominal_b) "P.Value" else "adj.P.Val"
    
    assays_a <- unique(df_a$assay); assays_b <- unique(df_b$assay)
    total_mat <- matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b))
    
    for(a in assays_a) {
        genes_a <- df_a$ID[df_a$assay == a & df_a[[p_col_a]] < 0.05]
        for(b in assays_b) {
            genes_b <- df_b$ID[df_b$assay == b & df_b[[p_col_b]] < 0.05]
            total_mat[a, b] <- length(intersect(genes_a, genes_b))
        }
    }
    
    ht <- Heatmap(total_mat, name = "Total Overlap", 
            border = TRUE,
            col = colorRamp2(c(0, max(total_mat)), c("white", "purple")),
            cluster_rows = TRUE, cluster_columns = TRUE,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(total_mat[i, j], x, y, gp = gpar(fontsize = 10))
            })
    options(repr.plot.width = plot.width, repr.plot.height = plot.height, repr.plot.res = plot.res)
    draw(ht)
}

#' Calculate directional overlaps
#' @param df_a, df_b DataFrames (must contain ID, assay, P.Value/adj.P.Val, logFC)
calc_directional_mats <- function(df_a, df_b,
                                  nominal_a = TRUE, nominal_b = TRUE) {
    
    # Standardize to data frames
    df_a <- as.data.frame(df_a)
    df_b <- as.data.frame(df_b)
    
    assays_a <- unique(df_a$assay)
    assays_b <- unique(df_b$assay)
    
    # Define p-value columns based on nominal flag
    p_col_a <- if(nominal_a) "P.Value" else "adj.P.Val"
    p_col_b <- if(nominal_b) "P.Value" else "adj.P.Val"
    
    # Initialize storage
    mats <- list(
        up_jaccard = matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b)),
        dn_jaccard = matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b)),
        up_counts  = matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b)),
        dn_counts  = matrix(0, length(assays_a), length(assays_b), dimnames = list(assays_a, assays_b))
    )
    
    for (a in assays_a) {
        # Subset df_a once per outer loop for efficiency
        sub_a <- df_a[df_a$assay == a, ]
        
        for (b in assays_b) {
            # Subset df_b once per inner loop
            sub_b <- df_b[df_b$assay == b, ]
            
            # --- UP-REGULATED ---
            genes_a_up <- sub_a$ID[sub_a[[p_col_a]] < 0.05 & sub_a$logFC > 0]
            genes_b_up <- sub_b$ID[sub_b[[p_col_b]] < 0.05 & sub_b$logFC > 0]
            
            inter_up <- length(intersect(genes_a_up, genes_b_up))
            uni_up   <- length(union(genes_a_up, genes_b_up))
            
            mats$up_jaccard[a, b] <- if (uni_up > 0) inter_up / uni_up else 0
            mats$up_counts[a, b]  <- inter_up
            
            # --- DOWN-REGULATED ---
            genes_a_dn <- sub_a$ID[sub_a[[p_col_a]] < 0.05 & sub_a$logFC < 0]
            genes_b_dn <- sub_b$ID[sub_b[[p_col_b]] < 0.05 & sub_b$logFC < 0]
            
            inter_dn <- length(intersect(genes_a_dn, genes_b_dn))
            uni_dn   <- length(union(genes_a_dn, genes_b_dn))
            
            mats$dn_jaccard[a, b] <- if (uni_dn > 0) inter_dn / uni_dn else 0
            mats$dn_counts[a, b]  <- inter_dn
        }
    }
    
    message(sprintf("[%s] Directional matrices calculated for %d x %d assays", 
                    Sys.time(), length(assays_a), length(assays_b)))
    
    return(mats)
}

COEF <- "Diff"
organoid_df <- topTable(readRDS('2026_03_03_Scanpy_020_2dpass_annotated_CLASS_cova_DiffPDvsCTRL_1.9.rds'), coef = COEF, number = Inf)
tissue_df <- topTable(readRDS('Lyra_data/dreamlet_class_upitt.rds'), coef = COEF, number = Inf)
#tissue_df <- readRDS('Lyra_data/subclass_upitt_dreamlet_toptable.rds')
tissue_sets <- get_gene_list(tissue_df, "T_")
organoid_sets <- get_gene_list(organoid_df, "O_")

tissue_df

# # Nominal - ADJusted

nominal_a = TRUE
nominal_b = TRUE
# 1. Prepare the matrices (using the logic from the previous step)
res_mats <- calc_directional_mats(organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b)
plt_total_overlap(organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b)

plt_total_vs_directional(res_mats,organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b,plot.width = 8)

plt_overalap_degs(res_mats,nominal_a = nominal_a,nominal_b = nominal_b,plot.width = 8)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "t"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.05         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s \n Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)
# pheatmap(
#   cluster_cor,
#   main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
#   color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
#   breaks = seq(-limit, limit, length.out = 101),
#   display_numbers = TRUE, number_format = "%.2f",
#   cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
# )
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize_number = 8,
  fontsize_main = 10
)
pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "t"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.05         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s | Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)

options(repr.plot.width = 6, repr.plot.height = 3, repr.plot.res = 200)
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f",
  cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
)

pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)

# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "union"      # "union" or "intersect"
CHOSEN_STAT   <- "logFC"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.01         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s | Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f",
  cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
)

pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "logFC"          # "t" or "logFC"
USE_CONC      <- FALSE        # TRUE or FALSE
PVAL_THRES    <- 0.01         # 0.05 or 0.01

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC,
                                    pval_thres = PVAL_THRES)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s | Genes: %d | PVAL_THRES: %.2f", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features),PVAL_THRES)

limit <- max(abs(cluster_cor), na.rm = TRUE)
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f",
  cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
)

pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "intersect"      # "union" or "intersect"
CHOSEN_STAT   <- "logFC"          # "t" or "logFC"
USE_CONC      <- FALSE         # TRUE or FALSE

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s | Genes: %d", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features))

limit <- max(abs(cluster_cor), na.rm = TRUE)
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f",
  cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
)

pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)


# ADJUST THESE TOGGLES:
CHOSEN_METHOD <- "union"      # "union" or "intersect"
CHOSEN_STAT   <- "t"          # "t" or "logFC"
USE_CONC      <- FALSE         # TRUE or FALSE

org_df_flat <- as.data.frame(organoid_df)
tiss_df_flat <- as.data.frame(tissue_df)

final_features <- get_feature_genes(org_df_flat, tiss_df_flat, 
                                    method = CHOSEN_METHOD, 
                                    stat_type = CHOSEN_STAT, 
                                    use_concordance = USE_CONC)

log_info(sprintf("Features selected: %d genes.", length(final_features)))

mat_org <- prep_cluster_matrix(org_df_flat, final_features, CHOSEN_STAT)
mat_tiss <- prep_cluster_matrix(tiss_df_flat, final_features, CHOSEN_STAT)

common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")
n_overlaps  <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# --- 3. DYNAMIC VISUALIZATION ---

# Dynamic Subtitle Generation
plot_subtitle <- sprintf("Selection: %s | Stat: %s | Concordance: %s | Genes: %d", 
                         toupper(CHOSEN_METHOD), CHOSEN_STAT, USE_CONC, length(final_features))

limit <- max(abs(cluster_cor), na.rm = TRUE)
pheatmap(
  cluster_cor,
  main = paste0("Cluster Spearman Correlation\n", plot_subtitle),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f",
  cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45
)

pheatmap(
  n_overlaps,
  main = paste0("Shared Gene Overlap Counts\n", plot_subtitle),
  display_numbers = TRUE, number_format = "%.0f",
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45
)



library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}


# 1. Background: Genes physically present in both dataframes
all_ids_org <- unique(organoid_df$ID)
all_ids_tiss <- unique(tissue_df$ID)
common_background <- intersect(all_ids_org, all_ids_tiss)

# 2. Genes of Interest: Union of significant genes (Nominal P < 0.05)
sig_ids_org <- organoid_df$ID[organoid_df$P.Value < 0.05]
sig_ids_tiss <- tissue_df$ID[tissue_df$P.Value < 0.05]
union_sig <- union(sig_ids_org, sig_ids_tiss)

# 3. Final Gene Set: Significant in at least one, but measured in both
final_features <- intersect(union_sig, common_background)


# 1. Prepare Matrices with NAs intact
prep_matrix_na <- function(df, target_genes) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_genes) %>%
    dplyr::group_by(assay, ID) %>%
    dplyr::summarize(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    # Not using values_fill = 0; letting missing combinations become NA
    tidyr::pivot_wider(names_from = assay, values_from = logFC) %>%
    tibble::column_to_rownames("ID")
}

log_info("Pivoting dataframes and retaining NA values for missing genes...")
mat_org <- prep_matrix_na(organoid_df, final_features)
mat_tiss <- prep_matrix_na(tissue_df, final_features)

# 2. Align Rows
# We ensure both matrices have the exact same genes, even if completely NA in one
all_genes <- sort(union(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- mat_org[all_genes, , drop = FALSE]
mat_tiss_final <- mat_tiss[all_genes, , drop = FALSE]

log_info("Calculating Spearman correlation using pairwise complete observations...")

# 3. Pairwise Correlation
# This calculates the correlation pair-by-pair, ignoring NAs only for the specific clusters being compared
cor_matrix <- cor(mat_org_final, mat_tiss_final, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

# 4. Zero-Centered Visualization
limit <- max(abs(cor_matrix), na.rm = TRUE)
palette_length <- 100
my_color <- colorRampPalette(c("#4575b4", "white", "#d73027"))(palette_length)
my_breaks <- seq(-limit, limit, length.out = palette_length + 1)

pheatmap(
  cor_matrix,
  color = my_color,
  breaks = my_breaks,
  main = "Spearman Correlation (Pairwise Complete Obs)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)

library(dplyr)
library(pheatmap)

#' Calculate the number of overlapping non-NA observations between two matrices
#' @param mat_x First matrix (genes as rows, clusters as columns)
#' @param mat_y Second matrix (genes as rows, clusters as columns)
#' @return A matrix of counts (clusters_x rows, clusters_y columns)
get_pairwise_n_matrix <- function(mat_x, mat_y) {
  
  cat(sprintf("[%s] INFO: Calculating pairwise overlapping gene counts...\n", Sys.time()), file = stderr())
  
  # 1. Create binary indicator matrices (TRUE = 1, FALSE = 0)
  # TRUE if the gene was measured in that cluster, FALSE if it is NA
  indicator_x <- !is.na(mat_x)
  indicator_y <- !is.na(mat_y)
  
  # 2. Matrix cross-multiplication
  # t(indicator_x) has dimensions (Clusters X by Genes)
  # indicator_y has dimensions (Genes by Clusters Y)
  # The dot product perfectly tallies the shared 1s (overlapping genes)
  n_matrix <- t(indicator_x) %*% indicator_y
  
  # Ensure the object is a standard matrix
  n_matrix <- as.matrix(n_matrix)
  
  cat(sprintf("[%s] INFO: Pairwise N matrix generated (%d x %d)\n", 
              Sys.time(), nrow(n_matrix), ncol(n_matrix)), file = stderr())
  
  return(n_matrix)
}

# Example Pipeline Integration:
# Assuming mat_org_final and mat_tiss_final are your matrices with NAs intact
n_overlaps <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# View the exact overlap between specific clusters:

# Optional: Visualize the overlap density
pheatmap(
  n_overlaps,
  main = "Number of Shared Genes Used for Correlation",
  display_numbers = TRUE,
  number_format = "%.0f",  # Integer formatting
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

length(final_features)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

# 1. Background: Genes physically present in both dataframes
all_ids_org <- unique(organoid_df$ID)
all_ids_tiss <- unique(tissue_df$ID)
common_background <- intersect(all_ids_org, all_ids_tiss)

# 2. Genes of Interest: Union of significant genes (Nominal P < 0.05)
sig_ids_org <- organoid_df$ID[organoid_df$P.Value < 0.05]
sig_ids_tiss <- tissue_df$ID[tissue_df$P.Value < 0.05]
union_sig <- union(sig_ids_org, sig_ids_tiss)

# Focus on genes that show signal in BOTH systems
shared_sig <- intersect(sig_ids_org, sig_ids_tiss)
final_features <- intersect(shared_sig, common_background)



# 1. Prepare Matrices with NAs intact
prep_matrix_na <- function(df, target_genes) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_genes) %>%
    dplyr::group_by(assay, ID) %>%
    dplyr::summarize(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    # Not using values_fill = 0; letting missing combinations become NA
    tidyr::pivot_wider(names_from = assay, values_from = logFC) %>%
    tibble::column_to_rownames("ID")
}

log_info("Pivoting dataframes and retaining NA values for missing genes...")
mat_org <- prep_matrix_na(organoid_df, final_features)
mat_tiss <- prep_matrix_na(tissue_df, final_features)

# 2. Align Rows
# We ensure both matrices have the exact same genes, even if completely NA in one
all_genes <- sort(union(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- mat_org[all_genes, , drop = FALSE]
mat_tiss_final <- mat_tiss[all_genes, , drop = FALSE]

log_info("Calculating Spearman correlation using pairwise complete observations...")

# 3. Pairwise Correlation
# This calculates the correlation pair-by-pair, ignoring NAs only for the specific clusters being compared
cor_matrix <- cor(mat_org_final, mat_tiss_final, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

# 4. Zero-Centered Visualization
limit <- max(abs(cor_matrix), na.rm = TRUE)
palette_length <- 100
my_color <- colorRampPalette(c("#4575b4", "white", "#d73027"))(palette_length)
my_breaks <- seq(-limit, limit, length.out = palette_length + 1)

pheatmap(
  cor_matrix,
  color = my_color,
  breaks = my_breaks,
  main = "Spearman Correlation (Pairwise Complete Obs)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)

length(final_features)

library(dplyr)
library(pheatmap)

#' Calculate the number of overlapping non-NA observations between two matrices
#' @param mat_x First matrix (genes as rows, clusters as columns)
#' @param mat_y Second matrix (genes as rows, clusters as columns)
#' @return A matrix of counts (clusters_x rows, clusters_y columns)
get_pairwise_n_matrix <- function(mat_x, mat_y) {
  
  cat(sprintf("[%s] INFO: Calculating pairwise overlapping gene counts...\n", Sys.time()), file = stderr())
  
  # 1. Create binary indicator matrices (TRUE = 1, FALSE = 0)
  # TRUE if the gene was measured in that cluster, FALSE if it is NA
  indicator_x <- !is.na(mat_x)
  indicator_y <- !is.na(mat_y)
  
  # 2. Matrix cross-multiplication
  # t(indicator_x) has dimensions (Clusters X by Genes)
  # indicator_y has dimensions (Genes by Clusters Y)
  # The dot product perfectly tallies the shared 1s (overlapping genes)
  n_matrix <- t(indicator_x) %*% indicator_y
  
  # Ensure the object is a standard matrix
  n_matrix <- as.matrix(n_matrix)
  
  cat(sprintf("[%s] INFO: Pairwise N matrix generated (%d x %d)\n", 
              Sys.time(), nrow(n_matrix), ncol(n_matrix)), file = stderr())
  
  return(n_matrix)
}

# Example Pipeline Integration:
# Assuming mat_org_final and mat_tiss_final are your matrices with NAs intact
n_overlaps <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# View the exact overlap between specific clusters:

# Optional: Visualize the overlap density
pheatmap(
  n_overlaps,
  main = "Number of Shared Genes Used for Correlation",
  display_numbers = TRUE,
  number_format = "%.0f",  # Integer formatting
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

organoid_df

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

# 1. Background and Significance Filtering
log_info("Filtering for shared significant genes...")

# Convert S4 DFrames to standard data.frames for dplyr compatibility
sig_org <- organoid_df %>% 
  as.data.frame() %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::select(ID, logFC_org = logFC)

sig_tiss <- tissue_df %>% 
  as.data.frame() %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::select(ID, logFC_tiss = logFC)

# Inner join for intersection of significance
shared_sig_df <- inner_join(sig_org, sig_tiss, by = "ID")

# 2. Directional Concordance Filter
# We keep genes where BOTH systems move in the same direction
concordant_genes_df <- shared_sig_df %>%
  dplyr::filter(sign(logFC_org) == sign(logFC_tiss))

concordant_ids <- concordant_genes_df$ID

log_info(sprintf("Shared Sig: %d | Concordant: %d", 
                 nrow(shared_sig_df), length(concordant_ids)))

# 3. Preparation Function (Handling S4 inputs)
prep_matrix_refined <- function(df, target_genes) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_genes) %>%
    dplyr::group_by(assay, ID) %>%
    dplyr::summarize(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = assay, values_from = logFC) %>%
    tibble::column_to_rownames("ID")
}

mat_org <- prep_matrix_refined(organoid_df, concordant_ids)
mat_tiss <- prep_matrix_refined(tissue_df, concordant_ids)

# 4. Alignment & Correlation
all_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- mat_org[all_genes, , drop = FALSE]
mat_tiss_final <- mat_tiss[all_genes, , drop = FALSE]

log_info("Calculating Spearman correlation on concordant shared signals...")
cor_matrix <- cor(mat_org_final, mat_tiss_final, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

# 5. Dynamic Heatmap Scaling
limit <- max(abs(cor_matrix), na.rm = TRUE)
if(is.na(limit) || limit == 0) limit <- 1 

my_color <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
my_breaks <- seq(-limit, limit, length.out = 101)

pheatmap(
  cor_matrix,
  color = my_color,
  breaks = my_breaks,
  main = "Spearman Correlation (Concordant Sig Genes)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)

library(dplyr)
library(pheatmap)

#' Calculate the number of overlapping non-NA observations between two matrices
#' @param mat_x First matrix (genes as rows, clusters as columns)
#' @param mat_y Second matrix (genes as rows, clusters as columns)
#' @return A matrix of counts (clusters_x rows, clusters_y columns)
get_pairwise_n_matrix <- function(mat_x, mat_y) {
  
  cat(sprintf("[%s] INFO: Calculating pairwise overlapping gene counts...\n", Sys.time()), file = stderr())
  
  # 1. Create binary indicator matrices (TRUE = 1, FALSE = 0)
  # TRUE if the gene was measured in that cluster, FALSE if it is NA
  indicator_x <- !is.na(mat_x)
  indicator_y <- !is.na(mat_y)
  
  # 2. Matrix cross-multiplication
  # t(indicator_x) has dimensions (Clusters X by Genes)
  # indicator_y has dimensions (Genes by Clusters Y)
  # The dot product perfectly tallies the shared 1s (overlapping genes)
  n_matrix <- t(indicator_x) %*% indicator_y
  
  # Ensure the object is a standard matrix
  n_matrix <- as.matrix(n_matrix)
  
  cat(sprintf("[%s] INFO: Pairwise N matrix generated (%d x %d)\n", 
              Sys.time(), nrow(n_matrix), ncol(n_matrix)), file = stderr())
  
  return(n_matrix)
}

# Example Pipeline Integration:
# Assuming mat_org_final and mat_tiss_final are your matrices with NAs intact
n_overlaps <- get_pairwise_n_matrix(mat_org_final, mat_tiss_final)

# View the exact overlap between specific clusters:

# Optional: Visualize the overlap density
pheatmap(
  n_overlaps,
  main = "Number of Shared Genes Used for Correlation",
  display_numbers = TRUE,
  number_format = "%.0f",  # Integer formatting
  color = colorRampPalette(c("white", "forestgreen"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

length(all_genes)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

# 1. Background and Significance Filtering (Intersection + Concordance)
log_info("Filtering for shared significant genes using t-statistic concordance...")

# Ensure standard data.frames
df_org_clean <- as.data.frame(organoid_df)
df_tiss_clean <- as.data.frame(tissue_df)

# Identify genes significant in BOTH (using P-value for selection)
sig_org <- df_org_clean %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::select(ID, t_org = t)

sig_tiss <- df_tiss_clean %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::select(ID, t_tiss = t)

# Inner join for intersection
shared_sig_df <- inner_join(sig_org, sig_tiss, by = "ID")

# Apply Concordance Filter based on t-statistic direction
concordant_ids <- shared_sig_df %>%
  dplyr::filter(sign(t_org) == sign(t_tiss)) %>%
  dplyr::pull(ID)

log_info(sprintf("Intersection Sig: %d | Concordant (t): %d", 
                 nrow(shared_sig_df), length(concordant_ids)))

# 2. Preparation Function (using t instead of logFC)
prep_matrix_t <- function(df, target_genes) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_genes) %>%
    dplyr::group_by(assay, ID) %>%
    # Use mean of t if there are technical replicates/multiple probes
    dplyr::summarize(t_stat = mean(t, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = assay, values_from = t_stat) %>%
    tibble::column_to_rownames("ID")
}

mat_org <- prep_matrix_t(df_org_clean, concordant_ids)
mat_tiss <- prep_matrix_t(df_tiss_clean, concordant_ids)

# 3. Alignment & Correlation
all_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- mat_org[all_genes, , drop = FALSE]
mat_tiss_final <- mat_tiss[all_genes, , drop = FALSE]

log_info("Calculating Spearman correlation on t-statistic ranks...")
cor_matrix <- cor(mat_org_final, mat_tiss_final, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

# 4. Visualization
limit <- max(abs(cor_matrix), na.rm = TRUE)
if(is.na(limit) || limit == 0) limit <- 1 

my_color <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
my_breaks <- seq(-limit, limit, length.out = 101)

pheatmap(
  cor_matrix,
  color = my_color,
  breaks = my_breaks,
  main = "Spearman Correlation (Concordant Sig t-stats)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

# 1. Background and Significance Filtering
log_info("Collapsing duplicate IDs and filtering for shared sig genes...")

# Standardize and unique-ify IDs
process_sig_set <- function(df, stat_name) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(P.Value < 0.05) %>%
    # COLLAPSE DUPLICATES HERE to fix many-to-many error
    dplyr::group_by(ID) %>%
    dplyr::summarize(t_stat = mean(t, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(!!stat_name := t_stat)
}

sig_org <- process_sig_set(organoid_df, "t_org")
sig_tiss <- process_sig_set(tissue_df, "t_tiss")

# Now inner_join will be 1-to-1
shared_sig_df <- inner_join(sig_org, sig_tiss, by = "ID")

# 2. Directional Concordance Filter (t-stat)
concordant_ids <- shared_sig_df %>%
  dplyr::filter(sign(t_org) == sign(t_tiss)) %>%
  dplyr::pull(ID)

log_info(sprintf("Intersection Sig: %d | Concordant (t): %d", 
                 nrow(shared_sig_df), length(concordant_ids)))

# 3. Preparation Function for Matrices
prep_matrix_t <- function(df, target_genes) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_genes) %>%
    dplyr::group_by(assay, ID) %>%
    dplyr::summarize(t_stat = mean(t, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = assay, values_from = t_stat) %>%
    tibble::column_to_rownames("ID")
}

mat_org <- prep_matrix_t(organoid_df, concordant_ids)
mat_tiss <- prep_matrix_t(tissue_df, concordant_ids)

# 4. Alignment & Correlation
all_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- mat_org[all_genes, , drop = FALSE]
mat_tiss_final <- mat_tiss[all_genes, , drop = FALSE]

log_info("Calculating Spearman correlation on unique t-statistic ranks...")
cor_matrix <- cor(mat_org_final, mat_tiss_final, 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

# 5. Visualization
limit <- max(abs(cor_matrix), na.rm = TRUE)
if(is.na(limit) || limit == 0) limit <- 1 

my_color <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)
my_breaks <- seq(-limit, limit, length.out = 101)

pheatmap(
  cor_matrix,
  color = my_color,
  breaks = my_breaks,
  main = "Spearman Correlation (Concordant Sig t-stats)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = 45
)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

# 1. Clean and Summarize t-statistics per Cluster
# We must ensure 1 row per Gene ID within each cluster (assay)
prepare_stat_matrix <- function(df, target_ids) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(ID %in% target_ids) %>%
    dplyr::group_by(assay, ID) %>%
    # Average t-stats if multiple entries exist for the same gene in one cluster
    dplyr::summarize(t_stat = mean(t, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = assay, values_from = t_stat) %>%
    tibble::column_to_rownames("ID")
}

# 2. Define the Gene Set (Intersection + Concordance)
# We calculate concordance based on the overall trend to pick our "features"
df_org_all <- as.data.frame(organoid_df)
df_tiss_all <- as.data.frame(tissue_df)

# Get shared significant genes
sig_ids_org <- unique(df_org_all$ID[df_org_all$P.Value < 0.05])
sig_ids_tiss <- unique(df_tiss_all$ID[df_tiss_all$P.Value < 0.05])
shared_sig <- intersect(sig_ids_org, sig_ids_tiss)

# Filter for directional agreement (using means across clusters to determine trend)
avg_t_org <- df_org_all %>% filter(ID %in% shared_sig) %>% group_by(ID) %>% summarize(m = mean(t))
avg_t_tiss <- df_tiss_all %>% filter(ID %in% shared_sig) %>% group_by(ID) %>% summarize(m = mean(t))

concordant_ids <- inner_join(avg_t_org, avg_t_tiss, by = "ID") %>%
  filter(sign(m.x) == sign(m.y)) %>%
  pull(ID)

log_info(sprintf("Using %d concordant genes for cluster correlation.", length(concordant_ids)))

# 3. Create the two matrices
mat_org <- prepare_stat_matrix(df_org_all, concordant_ids)
mat_tiss <- prepare_stat_matrix(df_tiss_all, concordant_ids)

# 4. Align Rows (Genes)
# Spearman needs the same genes in the same order in both matrices
common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

# 5. Correlate Clusters (Columns of mat_org vs Columns of mat_tiss)
log_info("Running Spearman correlation between Organoid clusters and Tissue clusters...")
cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")

# 6. Visualization
limit <- max(abs(cluster_cor), na.rm = TRUE)
if(is.na(limit) || limit == 0) limit <- 1

pheatmap(
  cluster_cor,
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  main = "Cluster Similarity: Organoid vs Tissue (t-stat Spearman)",
  display_numbers = TRUE,
  number_format = "%.2f",
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

log_info <- function(msg) {
  cat(sprintf("[%s] INFO: %s\n", Sys.time(), msg), file = stderr())
}

# 1. Standardize Inputs
df_org_all <- as.data.frame(organoid_df)
df_tiss_all <- as.data.frame(tissue_df)

# 2. Identify the Union of Significant Genes
sig_ids_org <- unique(df_org_all$ID[df_org_all$P.Value < 0.05])
sig_ids_tiss <- unique(df_tiss_all$ID[df_tiss_all$P.Value < 0.05])
union_sig_ids <- union(sig_ids_org, sig_ids_tiss)

# 3. Filter for Presence in Both and Directional Concordance
# Calculate global trend per gene to determine direction
trend_org <- df_org_all %>% 
  filter(ID %in% union_sig_ids) %>% 
  group_by(ID) %>% 
  summarize(t_org = mean(t, na.rm = TRUE), .groups = "drop")

trend_tiss <- df_tiss_all %>% 
  filter(ID %in% union_sig_ids) %>% 
  group_by(ID) %>% 
  summarize(t_tiss = mean(t, na.rm = TRUE), .groups = "drop")

# Join trends and keep only concordant genes present in both
concordant_features <- inner_join(trend_org, trend_tiss, by = "ID") %>%
  filter(sign(t_org) == sign(t_tiss)) %>%
  pull(ID)

log_info(sprintf("Final feature set size (Concordant Union): %d genes", length(concordant_features)))

# 4. Prepare Matrices (Genes in Rows, Clusters in Columns)
prepare_cluster_matrix <- function(df, target_ids) {
  df %>%
    as.data.frame() %>%
    filter(ID %in% target_ids) %>%
    group_by(assay, ID) %>%
    # Average t-stats to handle many-to-many duplicates within a cluster
    summarize(t_stat = mean(t, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = assay, values_from = t_stat) %>%
    column_to_rownames("ID")
}

mat_org <- prepare_cluster_matrix(df_org_all, concordant_features)
mat_tiss <- prepare_cluster_matrix(df_tiss_all, concordant_features)

# 5. Align Rows and Run Spearman Correlation
common_genes <- sort(intersect(rownames(mat_org), rownames(mat_tiss)))
mat_org_final <- as.matrix(mat_org[common_genes, , drop = FALSE])
mat_tiss_final <- as.matrix(mat_tiss[common_genes, , drop = FALSE])

log_info("Correlating clusters based on t-statistic fingerprints...")
cluster_cor <- cor(mat_org_final, mat_tiss_final, method = "spearman", use = "pairwise.complete.obs")

# 6. Visualization
limit <- max(abs(cluster_cor), na.rm = TRUE)
if(is.na(limit) || limit == 0) limit <- 1

pheatmap(
  cluster_cor,
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-limit, limit, length.out = 101),
  main = "Cluster Correlation: Concordant Union (t-stats)",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  angle_col = 45
)

any(is.na(tissue_df[tissue_df$ID %in% final_features, 'logFC']))

# Check if a specific gene exists in all clusters
# If this number is less than the total number of unique assays, 
# you will get NAs during the pivot.
gene_counts = length(organoid_df[organoid_df$ID == "A2M", 'assay'])
total_assays = length(unique(organoid_df$assay))
gene_counts
total_assays

organoid_df[organoid_df$ID == "A2M", ]

nominal_a = TRUE
nominal_b = TRUE
# 1. Prepare the matrices (using the logic from the previous step)
res_mats <- calc_directional_mats(organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b)

plt_total_overlap(organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b)

plt_total_vs_directional(res_mats,organoid_df,tissue_df,nominal_a = nominal_a,nominal_b = nominal_b,plot.width = 12.5,plot.height = 4)

plt_overalap_degs(res_mats,nominal_a = nominal_a,nominal_b = nominal_b)

# 2. Combine and build intersection matrix
combined_sets <- c(tissue_sets, organoid_sets)
m <- make_comb_mat(combined_sets)

# 3. Filter for significant overlaps to keep the plot clean
# We only show intersections with at least 10 genes
m_filtered <- m[comb_size(m) >= 5]

# 4. Plot
message(sprintf("[%s] INFO: Plotting UpSet for %d intersections", Sys.time(), ncol(m_filtered)))
options(repr.plot.width = 12, repr.plot.height = 10, repr.plot.res = 200)
UpSet(
  m_filtered,
  set_order = order(set_size(m_filtered), decreasing = TRUE),
  comb_order = order(comb_size(m_filtered), decreasing = TRUE),
  top_annotation = upset_top_annotation(m_filtered, add_numbers = TRUE),
  right_annotation = upset_right_annotation(m_filtered, add_numbers = TRUE),
  column_title = "Gene Set Intersections: Human Tissue vs Organoids"
)
