#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd gatk
require_cmd bcftools
require_cmd samtools
require_file "$TARGET_BED"
require_file "$ON_TARGET_TRUTH_VCF"

UNFILT="$VCF_DIR/mix124_vs_mix01_200M.unfiltered_TargetRegionOnly.vcf.gz"
FILT="$VCF_DIR/mix124_vs_mix01_200M.filtered_TargetRegionOnly.vcf.gz"

gatk --java-options "-Xmx${JAVA_XMX}" Mutect2 \
  -R "$REF_FASTA" \
  -I "$BAM_DIR/mix124_200M.markdup.bam" \
  -tumor mix124 \
  -I "$BAM_DIR/mix01_200M.markdup.bam" \
  -normal mix01 \
  --germline-resource "$GNOMAD_RESOURCE" \
  --panel-of-normals "$PON_RESOURCE" \
  -L "$TARGET_BED" \
  -O "$UNFILT"

[[ -f "${UNFILT}.tbi" ]] || gatk IndexFeatureFile -I "$UNFILT"
[[ -f "${UNFILT}.stats" ]] || { echo "Missing Mutect2 stats file: ${UNFILT}.stats" >&2; exit 1; }

gatk FilterMutectCalls \
  -R "$REF_FASTA" \
  -V "$UNFILT" \
  -O "$FILT"

echo -n "Tumor-normal PASS calls: "
bcftools view -f PASS -H "$FILT" | wc -l

bcftools isec -p "$VCF_DIR/truth_intersect_tumor_normal" \
  "$FILT" \
  "$ON_TARGET_TRUTH_VCF"

echo "Tumor-normal truth recovery:"
echo -n "Calls unique: "
bcftools view -H "$VCF_DIR/truth_intersect_tumor_normal/0000.vcf" | wc -l
echo -n "Truth unique (missed): "
bcftools view -H "$VCF_DIR/truth_intersect_tumor_normal/0001.vcf" | wc -l
echo -n "Shared true positives: "
bcftools view -H "$VCF_DIR/truth_intersect_tumor_normal/0002.vcf" | wc -l

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AF]\n' \
  "$VCF_DIR/truth_intersect_tumor_normal/0002.vcf" > "$VCF_DIR/recovered_truth_tumor_normal.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  "$VCF_DIR/truth_intersect_tumor_normal/0001.vcf" > "$VCF_DIR/missed_truth_tumor_normal.tsv"

awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $1":"$2":"$3">"$4}' \
  "$VCF_DIR/missed_truth_tumor_normal.tsv" > "$VCF_DIR/missed_truth_tumor_normal.bed"

samtools depth -b "$VCF_DIR/missed_truth_tumor_normal.bed" \
  "$BAM_DIR/mix124_200M.markdup.bam" > "$VCF_DIR/missed_truth_tumor_normal.depth.tsv"
