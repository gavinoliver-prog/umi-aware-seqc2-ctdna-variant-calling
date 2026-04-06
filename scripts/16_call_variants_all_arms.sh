#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Variant calling for all four pipeline arms.
#
# Arm A: Mutect2 tumor-only  + LoFreq tumor-only  (standard markdup BAM)
# Arm B: Mutect2 tumor-normal                      (standard markdup BAMs)
# Arm C: Mutect2 tumor-only  + LoFreq tumor-only  (consensus BAM)
# Arm D: Mutect2 tumor-normal                      (consensus BAMs)
#
# All calling restricted to $TARGET_BED with -L (Mutect2) / -l (LoFreq).
# gnomAD and PON used when available, skipped gracefully when absent.
#
# Output per arm directory:
#   mutect2.unfiltered.vcf.gz
#   mutect2.filtered.vcf.gz     (PASS variants ready for evaluation)
#   lofreq.vcf.gz               (Arms A and C only)
#
# Contamination tables written to the tumor BAM's directory so they are
# computed once and shared across arms that use the same tumor BAM:
#   ARM_A_DIR/tumor.{pileup,contamination}.table  — used by Arms A and B
#   ARM_C_DIR/tumor.{pileup,contamination}.table  — used by Arms C and D

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd gatk
require_cmd lofreq
require_cmd bcftools
require_file "$REF_FASTA"
require_file "$TARGET_BED"

# Build optional resource arguments for Mutect2.
# gnomAD germline resource improves germline variant filtering accuracy.
# PON suppresses recurrent sequencing artefacts seen across normal samples.
# Both improve specificity but are not required — Mutect2 runs without them.
MUTECT2_EXTRA=()
if [[ -f "$GNOMAD_RESOURCE" ]]; then
  MUTECT2_EXTRA+=(--germline-resource "$GNOMAD_RESOURCE")
  log "gnomAD resource found — will use for germline filtering and contamination."
else
  log "WARNING: gnomAD not found at $GNOMAD_RESOURCE — proceeding without it."
fi
if [[ -f "$PON_RESOURCE" ]]; then
  MUTECT2_EXTRA+=(--panel-of-normals "$PON_RESOURCE")
  log "PON found — will use for artefact filtering."
else
  log "WARNING: PON not found at $PON_RESOURCE — proceeding without it."
fi

# -------------------------------------------------------------------------
# run_contamination_estimate
#
# Runs GetPileupSummaries and CalculateContamination on the tumor BAM.
# Writes output to the BAM's own directory so it is shared across arms
# that use the same tumor BAM (Arms A+B share; Arms C+D share).
# Requires gnomAD resource — skips silently if unavailable.
# Idempotent: skips if contamination table already exists.
# -------------------------------------------------------------------------
run_contamination_estimate() {
  local tumor_bam="$1"
  local bam_dir
  bam_dir="$(dirname "$tumor_bam")"

  [[ -f "$GNOMAD_RESOURCE" ]] || {
    log "  No gnomAD resource — skipping contamination estimate."
    return 0
  }

  local pileup="$bam_dir/tumor.pileup.table"
  local contamination="$bam_dir/tumor.contamination.table"

  if [[ -f "$contamination" ]]; then
    log "  Contamination table already exists: $contamination"
    return 0
  fi

  log "  GetPileupSummaries..."
  gatk --java-options "-Xmx${JAVA_XMX}" GetPileupSummaries \
    -I "$tumor_bam" \
    -V "$GNOMAD_RESOURCE" \
    -L "$TARGET_BED" \
    -O "$pileup"

  log "  CalculateContamination..."
  gatk --java-options "-Xmx${JAVA_XMX}" CalculateContamination \
    -I "$pileup" \
    -O "$contamination"

  log "  Contamination table: $contamination"
}

# -------------------------------------------------------------------------
# mutect2_tumor_only: Mutect2 + contamination estimate + FilterMutectCalls
# -------------------------------------------------------------------------
mutect2_tumor_only() {
  local tumor_bam="$1"
  local out_dir="$2"

  local unfiltered="$out_dir/mutect2.unfiltered.vcf.gz"
  local filtered="$out_dir/mutect2.filtered.vcf.gz"

  if [[ -f "$filtered" ]]; then
    log "  $filtered exists — skipping."
    return 0
  fi

  require_file "$tumor_bam"

  log "  Mutect2 tumor-only: $(basename "$tumor_bam")..."
  gatk --java-options "-Xmx${JAVA_XMX}" Mutect2 \
    -R "$REF_FASTA" \
    -I "$tumor_bam" \
    --tumor-sample "tumor" \
    -L "$TARGET_BED" \
    "${MUTECT2_EXTRA[@]}" \
    -O "$unfiltered"

  # Contamination table lives in the BAM directory (shared across arms)
  run_contamination_estimate "$tumor_bam"

  local contam_args=()
  local bam_dir
  bam_dir="$(dirname "$tumor_bam")"
  [[ -f "$bam_dir/tumor.contamination.table" ]] && \
    contam_args+=(--contamination-table "$bam_dir/tumor.contamination.table")

  log "  FilterMutectCalls..."
  gatk --java-options "-Xmx${JAVA_XMX}" FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "$unfiltered" \
    "${contam_args[@]}" \
    -O "$filtered"

  log "  PASS calls: $(bcftools view -f PASS -H "$filtered" | wc -l)"
}

# -------------------------------------------------------------------------
# mutect2_tumor_normal: Mutect2 + contamination estimate + FilterMutectCalls
# -------------------------------------------------------------------------
mutect2_tumor_normal() {
  local tumor_bam="$1"
  local normal_bam="$2"
  local out_dir="$3"

  local unfiltered="$out_dir/mutect2.unfiltered.vcf.gz"
  local filtered="$out_dir/mutect2.filtered.vcf.gz"

  if [[ -f "$filtered" ]]; then
    log "  $filtered exists — skipping."
    return 0
  fi

  require_file "$tumor_bam"
  require_file "$normal_bam"

  log "  Mutect2 tumor-normal: $(basename "$tumor_bam") vs $(basename "$normal_bam")..."
  gatk --java-options "-Xmx${JAVA_XMX}" Mutect2 \
    -R "$REF_FASTA" \
    -I "$tumor_bam" \
    --tumor-sample "tumor" \
    -I "$normal_bam" \
    --normal-sample "normal" \
    -L "$TARGET_BED" \
    "${MUTECT2_EXTRA[@]}" \
    -O "$unfiltered"

  # Contamination table lives in the tumor BAM directory (may already exist
  # if tumor-only arm ran first with the same BAM)
  run_contamination_estimate "$tumor_bam"

  local contam_args=()
  local bam_dir
  bam_dir="$(dirname "$tumor_bam")"
  [[ -f "$bam_dir/tumor.contamination.table" ]] && \
    contam_args+=(--contamination-table "$bam_dir/tumor.contamination.table")

  log "  FilterMutectCalls..."
  gatk --java-options "-Xmx${JAVA_XMX}" FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "$unfiltered" \
    "${contam_args[@]}" \
    -O "$filtered"

  log "  PASS calls: $(bcftools view -f PASS -H "$filtered" | wc -l)"
}

# -------------------------------------------------------------------------
# lofreq_tumor_only: LoFreq SNV calling on a single tumor BAM
#
# LoFreq uses a Poisson-based model that directly accounts for base error
# rates, making it well-suited for low-VAF ctDNA calling as a complement
# to Mutect2's haplotype-aware approach.
#
# SNV calling only — indel calling requires lofreq indelqual to add
# per-read indel error qualities first, which would require writing a new
# BAM per arm. Omitted here; add indelqual step if indels are needed.
# -------------------------------------------------------------------------
lofreq_tumor_only() {
  local tumor_bam="$1"
  local out_dir="$2"

  local out_vcf="$out_dir/lofreq.vcf.gz"

  if [[ -f "$out_vcf" ]]; then
    log "  $out_vcf exists — skipping."
    return 0
  fi

  require_file "$tumor_bam"

  local tmp_vcf="$TMP_DIR/lofreq_$$.vcf"

  log "  LoFreq tumor-only: $(basename "$tumor_bam")..."
  lofreq call-parallel \
    --pp-threads "$THREADS" \
    -f "$REF_FASTA" \
    -l "$TARGET_BED" \
    --min-bq 20 \
    --min-alt-bq 20 \
    -o "$tmp_vcf" \
    "$tumor_bam"

  # Use -o flag (not shell redirection) to avoid bgzip stream corruption
  bcftools view -Oz -o "$out_vcf" "$tmp_vcf"
  bcftools index -t "$out_vcf"
  rm -f "$tmp_vcf"

  log "  LoFreq calls: $(bcftools view -H "$out_vcf" | wc -l)"
}

# =============================================================================
# Arm A: standard tumor-only
# =============================================================================
log "=== Arm A: standard tumor-only ==="
mutect2_tumor_only "$ARM_A_DIR/tumor.markdup.bam" "$ARM_A_DIR"
lofreq_tumor_only  "$ARM_A_DIR/tumor.markdup.bam" "$ARM_A_DIR"

# =============================================================================
# Arm B: standard tumor-normal
# Tumor BAM lives in ARM_A_DIR (shared with Arm A); normal BAM in ARM_B_DIR.
# The paired normal is Sample B — which is also the germline background in Ef.
# This shared germline signal causes over-suppression of real tumor variants
# (variants present in both tumor and normal get filtered as germline).
# Observing this empirically vs Arm A is a key result of this pipeline.
# Contamination table from ARM_A_DIR is reused if Arm A ran first.
# =============================================================================
log "=== Arm B: standard tumor-normal ==="
mutect2_tumor_normal \
  "$ARM_A_DIR/tumor.markdup.bam" \
  "$ARM_B_DIR/normal.markdup.bam" \
  "$ARM_B_DIR"

# =============================================================================
# Arm C: consensus tumor-only
# =============================================================================
log "=== Arm C: consensus tumor-only ==="
mutect2_tumor_only "$ARM_C_DIR/tumor.consensus.bam" "$ARM_C_DIR"
lofreq_tumor_only  "$ARM_C_DIR/tumor.consensus.bam" "$ARM_C_DIR"

# =============================================================================
# Arm D: consensus tumor-normal
# Tumor BAM is symlinked from ARM_C_DIR (created by script 14).
# Contamination table from ARM_C_DIR is reused if Arm C ran first.
# =============================================================================
log "=== Arm D: consensus tumor-normal ==="
mutect2_tumor_normal \
  "$ARM_D_DIR/tumor.consensus.bam" \
  "$ARM_D_DIR/normal.consensus.bam" \
  "$ARM_D_DIR"

log "Script 16 complete. Filtered VCFs:"
log "  Arm A: $ARM_A_DIR/mutect2.filtered.vcf.gz"
log "         $ARM_A_DIR/lofreq.vcf.gz"
log "  Arm B: $ARM_B_DIR/mutect2.filtered.vcf.gz"
log "  Arm C: $ARM_C_DIR/mutect2.filtered.vcf.gz"
log "         $ARM_C_DIR/lofreq.vcf.gz"
log "  Arm D: $ARM_D_DIR/mutect2.filtered.vcf.gz"
