#!/usr/bin/env bash
# =============================================================================
# 02_cellranger_annotate.sh
# Cell Ranger annotate — automated cell type annotation via 10x Genomics Cloud
#
# Author: Dimitrios Kyriakis
# Study: A53T alpha-synuclein iPSC-derived miniature brain organoid snRNA-seq
# Cell Ranger version: 9.0.1
# Annotation model: auto (GEX cell annotation, model selected automatically)
#
# Usage:
#   bash 02_cellranger_annotate.sh \
#       --count-dir  /path/to/01_cellranger_counts \
#       --outdir     /path/to/02_cellranger_annotate \
#       --token      /path/to/txg/credentials
#
# Required:
#   --count-dir  Directory produced by 01_cellranger_count.sh
#   --outdir     Output directory (will be created if absent)
#   --token      Path to 10x Genomics Cloud token file
#                (register at https://cloud.10xgenomics.com, then:
#                 cellranger cloud login --token-file /path/to/credentials)
#
# NOTE: This step requires internet access to 10x Genomics Cloud.
#       It is OPTIONAL — manual cell type annotation is performed in
#       03_qc_preprocessing_annotation.py (step 03).
#
# Output per sample:
#   <outdir>/<sample_id>/outs/cell_types.csv
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
COUNT_BASE=""
OUTPUT_BASE=""
TOKEN_PATH=""

usage() {
    sed -n '/^# Usage:/,/^# NOTE/{/^# NOTE/d; s/^# \{0,3\}//; p}' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --count-dir) COUNT_BASE="$2"; shift 2 ;;
        --outdir)    OUTPUT_BASE="$2"; shift 2 ;;
        --token)     TOKEN_PATH="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "[ERROR] Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$COUNT_BASE"   ]] && { echo "[ERROR] --count-dir is required"; usage; }
[[ -z "$OUTPUT_BASE"  ]] && { echo "[ERROR] --outdir is required";    usage; }
[[ -z "$TOKEN_PATH"   ]] && { echo "[ERROR] --token is required";     usage; }

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------
SAMPLES=(
    "Blanchard-wt-snRNAseq-a"
    "Blanchard-wt-snRNAseq-b"
    "Blanchard-A53T-E3-snRNAseq-a"
    "Blanchard-A53T-E3-snRNAseq-b"
    "Blanchard-A53T-E4-snRNAseq-a"
    "Blanchard-A53T-E4-snRNAseq-b"
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_BASE}"

# Required by Cell Ranger 9.x on some HPC systems
export MRO_DISK_SPACE_CHECK=disable

for SAMPLE_ID in "${SAMPLES[@]}"; do

    FLAG="${OUTPUT_BASE}/${SAMPLE_ID}_cellranger_annotate.done"

    if [[ -f "${FLAG}" ]]; then
        echo "[INFO] Skipping ${SAMPLE_ID} — already complete"
        continue
    fi

    MATRIX="${COUNT_BASE}/${SAMPLE_ID}/outs/filtered_feature_bc_matrix.h5"

    if [[ ! -f "${MATRIX}" ]]; then
        echo "[ERROR] Count matrix not found for ${SAMPLE_ID}: ${MATRIX}"
        echo "        Run 01_cellranger_count.sh first."
        exit 1
    fi

    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') | Running cellranger annotate: ${SAMPLE_ID}"

    cd "${OUTPUT_BASE}"

    cellranger annotate \
        --id="${SAMPLE_ID}" \
        --matrix="${MATRIX}" \
        --cell-annotation-model=auto \
        --tenx-cloud-token-path="${TOKEN_PATH}"

    # Remove large Cell Ranger intermediate files
    rm -rf "${OUTPUT_BASE}/${SAMPLE_ID}/CELLRANGER_ANNOTATE_CS"
    rm -rf "${OUTPUT_BASE}/${SAMPLE_ID}/_"*
    rm -rf "${OUTPUT_BASE}/${SAMPLE_ID}/extras"

    touch "${FLAG}"
    echo "[DONE] $(date '+%Y-%m-%d %H:%M:%S') | ${SAMPLE_ID}"

done

echo "[ALL DONE] Cell Ranger annotate complete for all ${#SAMPLES[@]} samples."
