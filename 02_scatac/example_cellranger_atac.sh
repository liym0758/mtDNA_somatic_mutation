#!/usr/bin/env bash
set -euo pipefail

# Example of the shared Cell Ranger ATAC count step.
: "${RUN_ID:?}"
: "${CUSTOM_REFERENCE:?}"
: "${FASTQ_DIR:?}"
: "${OUTPUT_DIR:?}"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
cellranger-atac count --id="$RUN_ID" --reference="$CUSTOM_REFERENCE" \
  --fastqs="$FASTQ_DIR" --localcores=16 --localmem=128
