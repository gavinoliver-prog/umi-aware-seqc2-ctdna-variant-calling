#!/usr/bin/env bash
set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$HOME/projects/ctdna2}"

export DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data}"
export FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/fastq}"
export SRA_DIR="${SRA_DIR:-$DATA_DIR/sra}"
export TMP_DIR="${TMP_DIR:-$DATA_DIR/tmp}"

export REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref}"
export TST170_DIR="${TST170_DIR:-$REF_DIR/IlluminaTST170_ref}"
export BRP_DIR="${BRP_DIR:-$REF_DIR/BRP2}"
export SAMPLEA_DIR="${SAMPLEA_DIR:-$REF_DIR/SampleA_ref}"

export RESULTS_DIR="${RESULTS_DIR:-$PROJECT_ROOT/results}"
export BAM_DIR="${BAM_DIR:-$RESULTS_DIR/bam}"
export QC_DIR="${QC_DIR:-$RESULTS_DIR/qc}"
export VCF_DIR="${VCF_DIR:-$RESULTS_DIR/vcf}"

export THREADS="${THREADS:-8}"
export JAVA_XMX="${JAVA_XMX:-24g}"

export REF_FASTA="${REF_FASTA:-$REF_DIR/Homo_sapiens_assembly38.fasta}"
export TARGET_BED="${TARGET_BED:-$BRP_DIR/BRP2_liftover_hg19tohg38.bed}"
export GNOMAD_RESOURCE="${GNOMAD_RESOURCE:-$REF_DIR/af-only-gnomad.hg38.vcf.gz}"
export PON_RESOURCE="${PON_RESOURCE:-$REF_DIR/1000g_pon.hg38.vcf.gz}"
export KNOWN_POSITIVES_VCF="${KNOWN_POSITIVES_VCF:-$SAMPLEA_DIR/KnownPositives_hg38.on_target.vcf.gz}"
export ON_TARGET_TRUTH_VCF="${ON_TARGET_TRUTH_VCF:-$SAMPLEA_DIR/KnownPositives_hg38.on_target.vcf.gz}"

export MIX01_RUN="${MIX01_RUN:-SRR13201012}"
export MIX124_RUN="${MIX124_RUN:-SRR13200999}"

# Semantic aliases for the 4-arm design
export TUMOR_RUN="${TUMOR_RUN:-SRR13200999}"   # Ef — 1:24 tumor:normal (~4% tumor fraction)
export NORMAL_RUN="${NORMAL_RUN:-SRR13201012}" # Bf — pure normal

# Four-arm output directories
export ARM_A_DIR="${ARM_A_DIR:-$RESULTS_DIR/arm_a_standard_tumor_only}"
export ARM_B_DIR="${ARM_B_DIR:-$RESULTS_DIR/arm_b_standard_tumor_normal}"
export ARM_C_DIR="${ARM_C_DIR:-$RESULTS_DIR/arm_c_consensus_tumor_only}"
export ARM_D_DIR="${ARM_D_DIR:-$RESULTS_DIR/arm_d_consensus_tumor_normal}"

# Truth set — known negatives for FP rate calculation
export KNOWN_NEGATIVES_BED="${KNOWN_NEGATIVES_BED:-$SAMPLEA_DIR/KnownNegatives_hg38.bed.gz}"

# Consensus calling parameters
export UMI_MIN_READS="${UMI_MIN_READS:-2}"
export UMI_MIN_BASE_QUALITY="${UMI_MIN_BASE_QUALITY:-20}"
export FGBIO_JAR="${FGBIO_JAR:-}"  # set by 11_install_missing_tools.sh

mkdir -p "$DATA_DIR" "$FASTQ_DIR" "$SRA_DIR" "$TMP_DIR" \
         "$REF_DIR" "$TST170_DIR" "$BRP_DIR" "$SAMPLEA_DIR" \
         "$RESULTS_DIR" "$BAM_DIR" "$QC_DIR" "$VCF_DIR" \
         "$ARM_A_DIR" "$ARM_B_DIR" "$ARM_C_DIR" "$ARM_D_DIR" \
         "$RESULTS_DIR/consensus_metrics" "$RESULTS_DIR/evaluation" \
         "$PROJECT_ROOT/scripts"

require_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || { echo "Missing command: $c" >&2; exit 1; }
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || { echo "Missing file: $f" >&2; exit 1; }
}
