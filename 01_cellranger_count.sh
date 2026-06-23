#!/usr/bin/env bash
# =============================================================================
# 01_cellranger_count.sh
# Cell Ranger count — FASTQ → filtered gene × barcode count matrices
#
# Author: Dimitrios Kyriakis
# Study: A53T alpha-synuclein iPSC-derived miniature brain organoid snRNA-seq
# Cell Ranger version: 9.0.1
# Reference genome: GRCh38-2024-A (10x Genomics pre-built human reference)
#
# Usage:
#   bash 01_cellranger_count.sh \
#       --reference  /path/to/refdata-gex-GRCh38-2024-A \
#       --fastqs     /path/to/fastq_files \
#       --outdir     /path/to/output/01_cellranger_counts \
#       --cores      8 \
#       --mem        64
#
# Required:
#   --reference  Path to Cell Ranger genome reference directory
#                Download: https://www.10xgenomics.com/support/software/cell-ranger/downloads
#   --fastqs     Directory containing raw FASTQ files (all 6 samples)
#   --outdir     Output directory (will be created if absent)
#
# Optional:
#   --cores      CPU cores per sample   [default: 8]
#   --mem        RAM in GB per sample   [default: 64]
#
# Samples processed (6 total — 2 replicates × 3 genotypes):
#   Blanchard-wt-snRNAseq-{a,b}         isogenic wild-type controls
#   Blanchard-A53T-E3-snRNAseq-{a,b}    A53T mutation, APOE-E3 background
#   Blanchard-A53T-E4-snRNAseq-{a,b}    A53T mutation, APOE-E4 background
#
# Output per sample:
#   <outdir>/<sample_id>/outs/filtered_feature_bc_matrix.h5
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
REFERENCE=""
FASTQ_PATH=""
OUTPUT_BASE=""
LOCALCORES=8
LOCALMEM=64

usage() {
    sed -n '/^# Usage:/,/^# =====/{/^# ====/d; s/^# \{0,3\}//; p}' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference) REFERENCE="$2"; shift 2 ;;
        --fastqs)    FASTQ_PATH="$2"; shift 2 ;;
        --outdir)    OUTPUT_BASE="$2"; shift 2 ;;
        --cores)     LOCALCORES="$2"; shift 2 ;;
        --mem)       LOCALMEM="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "[ERROR] Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$REFERENCE"  ]] && { echo "[ERROR] --reference is required"; usage; }
[[ -z "$FASTQ_PATH" ]] && { echo "[ERROR] --fastqs is required";    usage; }
[[ -z "$OUTPUT_BASE" ]] && { echo "[ERROR] --outdir is required";   usage; }

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------
SAMPLES=(
    "wt-snRNAseq-a"
    "wt-snRNAseq-b"
    "A53T-E3-snRNAseq-a"
    "A53T-E3-snRNAseq-b"
    "A53T-E4-snRNAseq-a"
    "A53T-E4-snRNAseq-b"
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_BASE}"

for SAMPLE_ID in "${SAMPLES[@]}"; do

    FLAG="${OUTPUT_BASE}/${SAMPLE_ID}_cellranger_count.done"

    if [[ -f "${FLAG}" ]]; then
        echo "[INFO] Skipping ${SAMPLE_ID} — already complete"
        continue
    fi

    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') | Running cellranger count: ${SAMPLE_ID}"

    cellranger count \
        --id="${SAMPLE_ID}" \
        --transcriptome="${REFERENCE}" \
        --fastqs="${FASTQ_PATH}" \
        --sample="${SAMPLE_ID}" \
        --localcores="${LOCALCORES}" \
        --localmem="${LOCALMEM}" \
        --create-bam=true \
        --output-dir="${OUTPUT_BASE}/${SAMPLE_ID}/"

    # Remove large Cell Ranger intermediate files to save disk space
    rm -rf "${OUTPUT_BASE}/${SAMPLE_ID}/SC_RNA_COUNTER_CS"
    rm -rf "${OUTPUT_BASE}/${SAMPLE_ID}/extras"

    touch "${FLAG}"
    echo "[DONE] $(date '+%Y-%m-%d %H:%M:%S') | ${SAMPLE_ID}"

done

echo "[ALL DONE] Cell Ranger count complete for all ${#SAMPLES[@]} samples."
