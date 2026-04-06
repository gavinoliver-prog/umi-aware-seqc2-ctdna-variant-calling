#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Collect and summarise consensus calling metrics.
#
# This script parses output files already written by scripts 13 and 14.
# No new alignments or fgbio runs — this should complete in seconds
# (the bcftools mpileup over ~53kb of known-negative space is the only
# computationally non-trivial step).
#
# Inputs (all pre-existing):
#   ARM_A_DIR/tumor.markdup_metrics.txt       — Picard MarkDuplicates metrics
#   ARM_A_DIR/tumor.flagstat.txt              — pre-consensus read counts
#   ARM_B_DIR/normal.markdup_metrics.txt
#   ARM_B_DIR/normal.flagstat.txt
#   ARM_C_DIR/tumor.grbu_family_sizes.txt     — fgbio GroupReadsByUmi histogram
#   ARM_C_DIR/tumor.flagstat.txt              — post-consensus read counts
#   ARM_D_DIR/normal.grbu_family_sizes.txt
#   ARM_D_DIR/normal.flagstat.txt
#   ARM_A_DIR/tumor.markdup.bam               — for pre-consensus error rate
#   ARM_C_DIR/tumor.consensus.bam             — for post-consensus error rate
#   ARM_B_DIR/normal.markdup.bam
#   ARM_D_DIR/normal.consensus.bam
#
# Outputs:
#   consensus_metrics/summary_metrics.tsv
#   consensus_metrics/family_size_distribution.tsv
#   consensus_metrics/error_rate_comparison.tsv

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd bcftools
require_cmd samtools
require_file "$REF_FASTA"
require_file "$KNOWN_NEGATIVES_BED"

METRICS_DIR="$RESULTS_DIR/consensus_metrics"
mkdir -p "$METRICS_DIR"

SUMMARY_TSV="$METRICS_DIR/summary_metrics.tsv"
FAMILY_TSV="$METRICS_DIR/family_size_distribution.tsv"
ERROR_TSV="$METRICS_DIR/error_rate_comparison.tsv"

if [[ -f "$SUMMARY_TSV" && -f "$FAMILY_TSV" && -f "$ERROR_TSV" ]]; then
  log "All output files exist — skipping."
  exit 0
fi

# -------------------------------------------------------------------------
# Helper: total reads from samtools flagstat
# Flagstat line: "N + 0 in total (QC-passed reads + QC-failed reads)"
# -------------------------------------------------------------------------
get_total_reads() {
  awk '/in total/ {print $1; exit}' "$1"
}

# -------------------------------------------------------------------------
# Helper: PERCENT_DUPLICATION from Picard MarkDuplicates metrics
# Locates column by header name — robust to column order changes across
# Picard/GATK versions.
# -------------------------------------------------------------------------
get_dup_rate() {
  awk '
    /^LIBRARY/ {
      for (i = 1; i <= NF; i++) if ($i == "PERCENT_DUPLICATION") col = i
      next
    }
    /^#/ || /^$/ { next }
    col { printf "%.4f\n", $col; exit }
  ' "$1"
}

# -------------------------------------------------------------------------
# Helper: error rate proxy via bcftools mpileup on known-negative regions
#
# Counts total and non-reference allele depths at base quality >= 20 across
# the 53,185 bp known-negative space. Non-ref rate is the proxy for the
# background error floor: pre-consensus this reflects raw sequencing error;
# post-consensus this reflects residual errors not corrected by family voting.
#
# Uses -T (--targets-file) rather than -R (--regions-file). Both restrict
# mpileup to specified regions, but -T reads the file sequentially and does
# NOT require a tabix index — important because KnownNegatives_hg38.bed.gz
# may be regular gzip rather than bgzip-compressed.
#
# Note on germline variants: Sample B germline variants present in both the
# Ef tumor and Bf normal samples will inflate the absolute error rate values.
# However, since both pre-consensus and post-consensus BAMs contain the same
# Sample B germline background, the PRE vs POST comparison remains valid —
# any reduction in the error rate proxy reflects consensus-calling error
# suppression, not germline variant removal. Absolute values should not be
# interpreted as a true sequencing error rate without further germline
# filtering.
# -------------------------------------------------------------------------
compute_error_rate() {
  local bam="$1"
  bcftools mpileup \
    -T "$KNOWN_NEGATIVES_BED" \
    -f "$REF_FASTA" \
    --min-BQ 20 \
    -a "AD" \
    -O u \
    "$bam" \
    | bcftools query -f '[%AD]\n' \
    | awk -F',' '
        NF >= 2 {
          ref += $1
          for (i = 2; i <= NF; i++) alt += $i
        }
        END {
          total = ref + alt
          if (total > 0) printf "%.8f\n", alt / total
          else           print "NA"
        }'
}

# -------------------------------------------------------------------------
# 1. Pre-consensus metrics — parse Picard MarkDuplicates and flagstat
# -------------------------------------------------------------------------
log "Parsing pre-consensus metrics..."

require_file "$ARM_A_DIR/tumor.markdup_metrics.txt"
require_file "$ARM_A_DIR/tumor.flagstat.txt"
require_file "$ARM_B_DIR/normal.markdup_metrics.txt"
require_file "$ARM_B_DIR/normal.flagstat.txt"

tumor_pre_reads=$(get_total_reads  "$ARM_A_DIR/tumor.flagstat.txt")
tumor_dup_rate=$(get_dup_rate      "$ARM_A_DIR/tumor.markdup_metrics.txt")
normal_pre_reads=$(get_total_reads "$ARM_B_DIR/normal.flagstat.txt")
normal_dup_rate=$(get_dup_rate     "$ARM_B_DIR/normal.markdup_metrics.txt")

log "  Tumor  pre-consensus reads: $tumor_pre_reads  dup rate: $tumor_dup_rate"
log "  Normal pre-consensus reads: $normal_pre_reads  dup rate: $normal_dup_rate"

# -------------------------------------------------------------------------
# 2. Post-consensus metrics — parse re-alignment flagstat from script 14
# -------------------------------------------------------------------------
log "Parsing post-consensus metrics..."

require_file "$ARM_C_DIR/tumor.flagstat.txt"
require_file "$ARM_D_DIR/normal.flagstat.txt"

tumor_post_reads=$(get_total_reads  "$ARM_C_DIR/tumor.flagstat.txt")
normal_post_reads=$(get_total_reads "$ARM_D_DIR/normal.flagstat.txt")

# Validate before arithmetic — a non-integer here means flagstat parsing failed
[[ "$tumor_post_reads"  =~ ^[0-9]+$ ]] \
  || { log "ERROR: tumor post-consensus read count is not a valid integer: '$tumor_post_reads'"; exit 1; }
[[ "$normal_post_reads" =~ ^[0-9]+$ ]] \
  || { log "ERROR: normal post-consensus read count is not a valid integer: '$normal_post_reads'"; exit 1; }

# Compression ratio: pre-consensus reads / post-consensus reads.
# Values >> 1 indicate aggressive duplicate collapse; ~1 indicates low duplication.
tumor_compression=$(awk  "BEGIN {printf \"%.2f\", $tumor_pre_reads  / $tumor_post_reads}")
normal_compression=$(awk "BEGIN {printf \"%.2f\", $normal_pre_reads / $normal_post_reads}")

log "  Tumor  post-consensus reads: $tumor_post_reads  compression: ${tumor_compression}x"
log "  Normal post-consensus reads: $normal_post_reads  compression: ${normal_compression}x"

# -------------------------------------------------------------------------
# 3. Family size distribution
# Parse the GroupReadsByUmi histogram written by script 14.
# fgbio format: family_size<tab>count<tab>fraction (with header line)
# -------------------------------------------------------------------------
log "Collating family size distributions..."

require_file "$ARM_C_DIR/tumor.grbu_family_sizes.txt"
require_file "$ARM_D_DIR/normal.grbu_family_sizes.txt"

{
  echo -e "sample\tfamily_size\tcount\tfraction"
  awk 'NR > 1 { print "tumor\t"  $0 }' "$ARM_C_DIR/tumor.grbu_family_sizes.txt"
  awk 'NR > 1 { print "normal\t" $0 }' "$ARM_D_DIR/normal.grbu_family_sizes.txt"
} > "$FAMILY_TSV"

log "Family size distribution: $FAMILY_TSV"

# -------------------------------------------------------------------------
# 4. Error rate proxy — bcftools mpileup on known negatives (~53kb, fast)
# -------------------------------------------------------------------------
log "Computing error rate proxy: tumor pre-consensus..."
tumor_pre_error=$(compute_error_rate "$ARM_A_DIR/tumor.markdup.bam")

log "Computing error rate proxy: tumor post-consensus..."
tumor_post_error=$(compute_error_rate "$ARM_C_DIR/tumor.consensus.bam")

log "Computing error rate proxy: normal pre-consensus..."
normal_pre_error=$(compute_error_rate "$ARM_B_DIR/normal.markdup.bam")

log "Computing error rate proxy: normal post-consensus..."
normal_post_error=$(compute_error_rate "$ARM_D_DIR/normal.consensus.bam")

{
  printf "sample\tcondition\terror_rate_proxy\n"
  printf "tumor\tpre_consensus\t%s\n"   "$tumor_pre_error"
  printf "tumor\tpost_consensus\t%s\n"  "$tumor_post_error"
  printf "normal\tpre_consensus\t%s\n"  "$normal_pre_error"
  printf "normal\tpost_consensus\t%s\n" "$normal_post_error"
} > "$ERROR_TSV"

log "Error rate comparison: $ERROR_TSV"

# -------------------------------------------------------------------------
# 5. Summary table
# -------------------------------------------------------------------------
{
  printf "metric\ttumor\tnormal\n"
  printf "pre_consensus_total_reads\t%s\t%s\n"  "$tumor_pre_reads"     "$normal_pre_reads"
  printf "markdup_duplicate_rate\t%s\t%s\n"     "$tumor_dup_rate"      "$normal_dup_rate"
  printf "post_consensus_total_reads\t%s\t%s\n" "$tumor_post_reads"    "$normal_post_reads"
  printf "compression_ratio\t%s\t%s\n"          "$tumor_compression"   "$normal_compression"
  printf "error_rate_pre_consensus\t%s\t%s\n"   "$tumor_pre_error"     "$normal_pre_error"
  printf "error_rate_post_consensus\t%s\t%s\n"  "$tumor_post_error"    "$normal_post_error"
} > "$SUMMARY_TSV"

log "Summary metrics: $SUMMARY_TSV"
log "Script 15 complete."
