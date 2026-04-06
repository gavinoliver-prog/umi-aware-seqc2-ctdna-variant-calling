#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd bwa
require_cmd samtools
require_cmd picard

sample_name() {
  case "$1" in
    "$MIX01_RUN") echo "mix01" ;;
    "$MIX124_RUN") echo "mix124" ;;
    *) echo "unknown" ;;
  esac
}

for depth in 5M 20M 100M 200M; do
  for run in "$MIX01_RUN" "$MIX124_RUN"; do
    sample="$(sample_name "$run")"
    prefix="${sample}_${depth}"
    r1="$FASTQ_DIR/sub_${depth}_${run}_R1.fastq.gz"
    r2="$FASTQ_DIR/sub_${depth}_${run}_R2.fastq.gz"

    require_file "$r1"
    require_file "$r2"

    bwa mem -t "$THREADS" \
      -R "@RG\tID:${run}_${depth}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1" \
      "$REF_FASTA" "$r1" "$r2" \
      | samtools sort -@ "$THREADS" -o "$BAM_DIR/${prefix}.coord.bam"

    samtools index "$BAM_DIR/${prefix}.coord.bam"

    samtools sort -n -@ "$THREADS" -o "$BAM_DIR/${prefix}.name.bam" "$BAM_DIR/${prefix}.coord.bam"
    samtools fixmate -m "$BAM_DIR/${prefix}.name.bam" "$BAM_DIR/${prefix}.fixmate.bam"
    samtools sort -@ "$THREADS" -o "$BAM_DIR/${prefix}.fixmate.coord.bam" "$BAM_DIR/${prefix}.fixmate.bam"
    samtools markdup -@ "$THREADS" "$BAM_DIR/${prefix}.fixmate.coord.bam" "$BAM_DIR/${prefix}.markdup.bam"
    samtools index "$BAM_DIR/${prefix}.markdup.bam"

    samtools flagstat "$BAM_DIR/${prefix}.markdup.bam" > "$QC_DIR/${prefix}.flagstat.txt"
    samtools stats "$BAM_DIR/${prefix}.markdup.bam" > "$QC_DIR/${prefix}.stats.txt"

    picard CollectInsertSizeMetrics \
      I="$BAM_DIR/${prefix}.markdup.bam" \
      O="$QC_DIR/${prefix}.insert_size_metrics.txt" \
      H="$QC_DIR/${prefix}.insert_size_histogram.pdf" \
      M=0.5
  done
done
