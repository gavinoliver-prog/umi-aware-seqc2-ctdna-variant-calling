#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd bcftools
require_cmd mosdepth
require_cmd samtools
require_file "$TARGET_BED"
require_file "$KNOWN_POSITIVES_VCF"

# Ensure truth VCF is indexed
if [ ! -f "${KNOWN_POSITIVES_VCF}.tbi" ]; then
  echo "Indexing truth VCF..."
  bcftools index -t "$KNOWN_POSITIVES_VCF"
fi

# Subset truth set to target regions
echo "Subsetting truth VCF to target regions..."
bcftools view \
  -R "$TARGET_BED" \
  "$KNOWN_POSITIVES_VCF" \
  -Oz \
  -o "$ON_TARGET_TRUTH_VCF"

bcftools index -f -t "$ON_TARGET_TRUTH_VCF"

# Coverage summary across depths
coverage_summary="$VCF_DIR/mix124.coverage_threshold_summary.tsv"
echo -e "sample\tgte_100x\tgte_250x\tgte_500x\tgte_1000x" > "$coverage_summary"

for depth in 20M 100M 200M; do
  prefix="$VCF_DIR/mix124_${depth}"

  mosdepth --by "$TARGET_BED" --thresholds 100,250,500,1000 \
    "$prefix" \
    "$BAM_DIR/mix124_${depth}.markdup.bam"

  zcat "${prefix}.thresholds.bed.gz" | awk -v sample="mix124_${depth}" '{
    len = $3 - $2
    sum100 += $5
    sum250 += $6
    sum500 += $7
    sum1000 += $8
    total += len
  } END {
    print sample "\t" sum100/total "\t" sum250/total "\t" sum500/total "\t" sum1000/total
  }' >> "$coverage_summary"
done

# Truth recovery summaries across depths
truth_summary_all="$VCF_DIR/mix124.truth_recovery_all_summary.tsv"
truth_summary_pass="$VCF_DIR/mix124.truth_recovery_pass_summary.tsv"

echo -e "sample\tcalls_unique\ttruth_missed\ttruth_recovered" > "$truth_summary_all"
echo -e "sample\tcalls_unique\ttruth_missed\ttruth_recovered" > "$truth_summary_pass"

for depth in 5M 20M 100M 200M; do
  filtered_vcf="$VCF_DIR/mix124_${depth}.filtered_TargetRegionOnly.vcf.gz"
  pass_vcf="$VCF_DIR/mix124_${depth}.filtered_PASS_TargetRegionOnly.vcf.gz"

  require_file "$filtered_vcf"

  # PASS-only subset
  bcftools view -f PASS -Oz \
    -o "$pass_vcf" \
    "$filtered_vcf"

  bcftools index -f -t "$pass_vcf"

  # ALL records from filtered VCF
  outdir_all="$VCF_DIR/truth_intersect_${depth}_all"
  bcftools isec -p "$outdir_all" \
    "$filtered_vcf" \
    "$ON_TARGET_TRUTH_VCF"

  calls_unique_all=$(bcftools view -H "$outdir_all/0000.vcf" | wc -l)
  truth_missed_all=$(bcftools view -H "$outdir_all/0001.vcf" | wc -l)
  truth_recovered_all=$(bcftools view -H "$outdir_all/0002.vcf" | wc -l)

  echo -e "mix124_${depth}\t${calls_unique_all}\t${truth_missed_all}\t${truth_recovered_all}" >> "$truth_summary_all"

  # PASS-only records
  outdir_pass="$VCF_DIR/truth_intersect_${depth}_pass"
  bcftools isec -p "$outdir_pass" \
    "$pass_vcf" \
    "$ON_TARGET_TRUTH_VCF"

  calls_unique_pass=$(bcftools view -H "$outdir_pass/0000.vcf" | wc -l)
  truth_missed_pass=$(bcftools view -H "$outdir_pass/0001.vcf" | wc -l)
  truth_recovered_pass=$(bcftools view -H "$outdir_pass/0002.vcf" | wc -l)

  echo -e "mix124_${depth}\t${calls_unique_pass}\t${truth_missed_pass}\t${truth_recovered_pass}" >> "$truth_summary_pass"
done

# Detailed recovered/missed truth analysis for 200M PASS only
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AF]\n' \
  "$VCF_DIR/truth_intersect_200M_pass/0002.vcf" > "$VCF_DIR/recovered_truth.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  "$VCF_DIR/truth_intersect_200M_pass/0001.vcf" > "$VCF_DIR/missed_truth.tsv"

awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $1":"$2":"$3">"$4}' \
  "$VCF_DIR/missed_truth.tsv" > "$VCF_DIR/missed_truth.bed"

samtools depth -b "$VCF_DIR/missed_truth.bed" \
  "$BAM_DIR/mix124_200M.markdup.bam" > "$VCF_DIR/missed_truth.depth.tsv"

awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $1":"$2":"$3">"$4}' \
  "$VCF_DIR/recovered_truth.tsv" > "$VCF_DIR/recovered_truth.bed"

samtools depth -b "$VCF_DIR/recovered_truth.bed" \
  "$BAM_DIR/mix124_200M.markdup.bam" > "$VCF_DIR/recovered_truth.depth.tsv"
