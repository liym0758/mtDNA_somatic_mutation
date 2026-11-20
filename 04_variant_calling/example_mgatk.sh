#!/usr/bin/env bash
set -euo pipefail

# Example of the shared species-level mgatk invocation.
: "${BAM_FILE:?}"
: "${BARCODE_FILE:?}"
: "${MT_FASTA:?}"
: "${SAMPLE_NAME:?}"
: "${OUTPUT_DIR:?}"
mkdir -p "$OUTPUT_DIR"
mgatk tenx -i "$BAM_FILE" -b "$BARCODE_FILE" -g "$MT_FASTA" \
  -bt CB -n "$SAMPLE_NAME" -o "$OUTPUT_DIR" --NMmax 2 \
  -q 30 -aq 30 -mb 1000 --proper-pairs --keep-qc-bams \
  --emit-base-qualities -c 12
