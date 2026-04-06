#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
source "$DIR/00_config.sh"

echo "=== ctDNA variant calling pipeline: full run ==="
echo "Project root: $PROJECT_ROOT"
echo ""

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
if [[ ! -f "$TARGET_BED" ]]; then
  echo "ERROR: TARGET_BED not found at $TARGET_BED" >&2
  echo "Expected: BRP2 panel BED at $BRP_DIR/BRP2_liftover_hg19tohg38.bed" >&2
  exit 1
fi

if [[ ! -f "$KNOWN_POSITIVES_VCF" ]]; then
  echo "ERROR: KNOWN_POSITIVES_VCF not found at $KNOWN_POSITIVES_VCF" >&2
  echo "Expected: SEQC2 BRP ground truth VCF in $SAMPLEA_DIR" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Phase 1: Environment setup and data preparation
# ---------------------------------------------------------------------------
echo ""
echo "--- Phase 1: Setup ---"

# Environment creation is handled by 01_setup_env.sh using environment.yml.
# 11_install_missing_tools.sh is superseded by this approach and kept only
# for reference (manual install into an existing environment).
bash "$DIR/01_setup_env.sh"
bash "$DIR/02_download_data.sh"

# Script 03 (subsample_fastq.sh) is optional — uncomment to run on a subset
# of reads for faster iteration during development.
# bash "$DIR/03_subsample_fastq.sh"

bash "$DIR/04_prepare_reference.sh"
bash "$DIR/05_fastqc_multiqc.sh"

# Scripts 06-10 are superseded by the 4-arm pipeline (scripts 12-18).
#   06_align_and_qc.sh             — replaced by 13_align_arms_ab.sh + 14_align_arms_cd.sh
#   07_call_mutect_tumor_only.sh   — replaced by 16_call_variants_all_arms.sh (Arm A)
#   08_coverage_and_truth.sh       — replaced by 17_evaluate_results.sh
#   09_call_mutect_tumor_normal.sh — replaced by 16_call_variants_all_arms.sh (Arm B)
#   10_compare_results.sh          — replaced by 17_evaluate_results.sh + 18_generate_report.sh

# ---------------------------------------------------------------------------
# Phase 2: Four-arm pipeline
# ---------------------------------------------------------------------------
echo ""
echo "--- Phase 2: Pipeline ---"

# Superseded by 01_setup_env.sh (environment.yml handles tool installation).
# Kept for reference — use to manually install fgbio/lofreq into an existing env.
# bash "$DIR/11_install_missing_tools.sh"

# Build synthetic UMI FASTQs (UMI-prefixed for Arms C/D; standard for Arms A/B)
bash "$DIR/12_group_fragments_by_coordinate.sh"

# Alignment for all arms.
# NOTE: Scripts 13 and 14 are independent and can be run in parallel:
#   bash "$DIR/13_align_arms_ab.sh" & bash "$DIR/14_align_arms_cd.sh" & wait
# Run sequentially here to keep log output readable in a single terminal.
bash "$DIR/13_align_arms_ab.sh"
bash "$DIR/14_align_arms_cd.sh"

# Quantify consensus calling effect (compression ratio, family size distribution)
bash "$DIR/15_collect_consensus_metrics.sh"

# Variant calling: Mutect2 + LoFreq across all four arms
bash "$DIR/16_call_variants_all_arms.sh"

# Evaluate sensitivity and FP rate against SEQC2 ground truth
bash "$DIR/17_evaluate_results.sh"

# Generate final report
bash "$DIR/18_generate_report.sh"

# ---------------------------------------------------------------------------
echo ""
echo "=== Pipeline complete ==="
echo "Key outputs:"
echo "  QC:          $RESULTS_DIR/qc/multiqc_report.html"
echo "  Consensus:   $RESULTS_DIR/consensus_metrics/summary_metrics.tsv"
echo "  Evaluation:  $RESULTS_DIR/evaluation/per_arm_summary.tsv"
echo "  Report:      $RESULTS_DIR/REPORT.md"
