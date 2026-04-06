#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd fastqc
require_cmd multiqc

mkdir -p "$QC_DIR/fastqc" "$QC_DIR/multiqc"

for run in "$MIX01_RUN" "$MIX124_RUN"; do
  fastqc -t "$THREADS" -o "$QC_DIR/fastqc" \
    "$FASTQ_DIR/sub_5M_${run}_R1.fastq.gz" \
    "$FASTQ_DIR/sub_5M_${run}_R2.fastq.gz"
done

multiqc "$QC_DIR/fastqc" -o "$QC_DIR/multiqc"
