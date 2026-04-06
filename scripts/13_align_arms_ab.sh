#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Align standard (non-consensus) FASTQs for Arms A and B.
#
# Arms A and B share the same tumor and normal BAMs — they differ only in
# how variant calling is run (tumor-only vs tumor-normal in script 16).
# The tumor BAM is written to ARM_A_DIR and referenced from there by both arms.
# The normal BAM is written to ARM_B_DIR (only needed for Arm B calling).
#
# Input:  on-target FASTQs produced by script 12 (${run}_ontarget_R1/R2.fastq.gz)
# Output: ARM_A_DIR/tumor.markdup.bam
#         ARM_B_DIR/normal.markdup.bam
#
# Duplicate marking uses Picard MarkDuplicates (via GATK4), which groups reads
# by fragment coordinates — the same grouping logic as script 12. This is the
# baseline against which consensus calling (Arms C/D) is compared.

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd bwa
require_cmd samtools
require_cmd gatk
require_file "$REF_FASTA"

align_sample() {
  local run="$1"
  local label="$2"    # sample name embedded in BAM read group; must match Mutect2 -tumor/-normal
  local out_bam="$3"  # full path to output markdup BAM

  local r1="$FASTQ_DIR/${run}_ontarget_R1.fastq.gz"
  local r2="$FASTQ_DIR/${run}_ontarget_R2.fastq.gz"
  require_file "$r1"
  require_file "$r2"

  local out_dir
  out_dir="$(dirname "$out_bam")"

  if [[ -f "${out_bam}.bai" ]]; then
    log "[$label] $out_bam already indexed — skipping."
    return 0
  fi

  local tmp_sorted="$TMP_DIR/${run}_ab_coord.bam"

  # BWA-MEM alignment.
  # -K 100000000: process reads in chunks of 100M bases for deterministic output.
  # Read group fields required by GATK: ID, SM, PL, LB.
  # SM must match the sample name passed to Mutect2 -tumor/-normal in script 16.
  log "[$label] Aligning to hg38..."
  bwa mem -t "$THREADS" -K 100000000 \
    -R "@RG\tID:${run}\tSM:${label}\tPL:ILLUMINA\tLB:lib1" \
    "$REF_FASTA" "$r1" "$r2" \
    | samtools sort -@ "$THREADS" -T "$TMP_DIR/${run}_ab_sort" -o "$tmp_sorted"

  # Picard MarkDuplicates: identifies duplicate fragments by coordinate,
  # flags them in the BAM but retains all reads (Mutect2 ignores flagged
  # duplicates internally). This is the standard-of-care baseline for
  # ctDNA — Arms A and B show what is achievable without consensus calling.
  #
  # OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500: NovaSeq uses a patterned flowcell
  # (ExAmp chemistry). Optical duplicates on patterned flowcells cluster at
  # greater pixel distances than on non-patterned flowcells (default 100).
  # 2500 is the GATK best-practices recommendation for NovaSeq.
  log "[$label] Running Picard MarkDuplicates..."
  gatk --java-options "-Xmx${JAVA_XMX}" MarkDuplicates \
    -I "$tmp_sorted" \
    -O "$out_bam" \
    -M "${out_dir}/${label}.markdup_metrics.txt" \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --TMP_DIR "$TMP_DIR"

  samtools index "$out_bam"
  rm -f "$tmp_sorted"

  # Alignment QC metrics
  log "[$label] Collecting metrics..."
  samtools flagstat -@ "$THREADS" "$out_bam" \
    > "${out_dir}/${label}.flagstat.txt"

  gatk --java-options "-Xmx${JAVA_XMX}" CollectInsertSizeMetrics \
    -I "$out_bam" \
    -O "${out_dir}/${label}.insert_size_metrics.txt" \
    -H "${out_dir}/${label}.insert_size_histogram.pdf"

  log "[$label] Done: $out_bam"
}

align_sample "$TUMOR_RUN"  "tumor"  "$ARM_A_DIR/tumor.markdup.bam"
align_sample "$NORMAL_RUN" "normal" "$ARM_B_DIR/normal.markdup.bam"

log "Script 13 complete."
log "  Tumor BAM  (Arms A+B): $ARM_A_DIR/tumor.markdup.bam"
log "  Normal BAM (Arm B):    $ARM_B_DIR/normal.markdup.bam"
