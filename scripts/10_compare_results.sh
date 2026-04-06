#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd bcftools

echo "Tumor-only PASS counts:"
for f in "$VCF_DIR"/mix01_*filtered_TargetRegionOnly.vcf.gz "$VCF_DIR"/mix124_*filtered_TargetRegionOnly.vcf.gz; do
  [[ -f "$f" ]] || continue
  printf "%s\t" "$(basename "$f")"
  bcftools view -f PASS -H "$f" | wc -l
done

echo
if [[ -f "$VCF_DIR/mix124_vs_mix01_200M.filtered_TargetRegionOnly.vcf.gz" ]]; then
  echo -n "Tumor-normal PASS count: "
  bcftools view -f PASS -H "$VCF_DIR/mix124_vs_mix01_200M.filtered_TargetRegionOnly.vcf.gz" | wc -l
fi

echo
if [[ -d "$VCF_DIR/truth_intersect" ]]; then
  echo "Tumor-only truth recovery (200M):"
  echo -n "Recovered: "
  bcftools view -H "$VCF_DIR/truth_intersect/0002.vcf" | wc -l
  echo -n "Missed: "
  bcftools view -H "$VCF_DIR/truth_intersect/0001.vcf" | wc -l
fi

if [[ -d "$VCF_DIR/truth_intersect_tumor_normal" ]]; then
  echo "Tumor-normal truth recovery (200M):"
  echo -n "Recovered: "
  bcftools view -H "$VCF_DIR/truth_intersect_tumor_normal/0002.vcf" | wc -l
  echo -n "Missed: "
  bcftools view -H "$VCF_DIR/truth_intersect_tumor_normal/0001.vcf" | wc -l
fi
