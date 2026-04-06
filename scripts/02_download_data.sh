#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd prefetch
require_cmd fasterq-dump
require_cmd pigz

download_and_extract() {
  local run="$1"
  prefetch "$run" --output-directory "$SRA_DIR" --max-size 50G

  fasterq-dump "$SRA_DIR/$run" \
    --split-files \
    --threads "$THREADS" \
    --outdir "$FASTQ_DIR" \
    --temp "$TMP_DIR"

  pigz -p "$THREADS" "$FASTQ_DIR/${run}_1.fastq"
  pigz -p "$THREADS" "$FASTQ_DIR/${run}_2.fastq"
}

#download_and_extract "$MIX01_RUN"
download_and_extract "$MIX124_RUN"
