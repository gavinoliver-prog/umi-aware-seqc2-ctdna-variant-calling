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
# Truth set: 228 known positives reduced to 86 somatic candidates
# Filtering: variants present in gnomAD at AF >= 0.001 excluded
# Rationale: 124 (54%) of original 228 are common germline SNPs
#            correctly filtered by Mutect2 with gnomAD resource
# Remaining 86 are not in gnomAD or extremely rare (AF < 0.001)
# All 86 are Tier 3/4 (<0.8% expected VAF in Ef, <4% in Df)
export KNOWN_POSITIVES_SOMATIC_VCF="${KNOWN_POSITIVES_SOMATIC_VCF:-$SAMPLEA_DIR/KnownPositives_hg38.somatic_only.vcf.gz}"
export KNOWN_POSITIVES_VCF="${KNOWN_POSITIVES_VCF:-$SAMPLEA_DIR/KnownPositives_hg38.somatic_only.vcf.gz}"
export ON_TARGET_TRUTH_VCF="${ON_TARGET_TRUTH_VCF:-$SAMPLEA_DIR/KnownPositives_hg38.somatic_only.vcf.gz}"

export MIX01_RUN="${MIX01_RUN:-SRR13201012}"
export MIX124_RUN="${MIX124_RUN:-SRR13200966}"

# Semantic aliases for the 4-arm design
# Df = Sample A diluted 1:4 into Sample B (SRR13200966)
# Expected VAF in Df = Sample A VAF / 5
# At 20% tumor fraction, Tier 1 variants appear at ~10-20% VAF
# This is within reliable detection range for standard calling
export TUMOR_RUN="${TUMOR_RUN:-SRR13200966}"   # Df — 1:4 tumor:normal (~20% tumor fraction)
export NORMAL_RUN="${NORMAL_RUN:-SRR13201012}" # Bf — pure normal (unchanged)
export DILUTION_FACTOR="${DILUTION_FACTOR:-5}" # 1:4 dilution = 1/5

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
export FGBIO_JAR="${FGBIO_JAR:-/home/ubuntu/tools/fgbio/fgbio-2.3.0.jar}"  # set by 11_install_missing_tools.sh

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
