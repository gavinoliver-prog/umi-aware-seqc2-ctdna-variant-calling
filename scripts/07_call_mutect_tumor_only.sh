#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd gatk
require_file "$TARGET_BED"
require_file "$GNOMAD_RESOURCE"
require_file "$PON_RESOURCE"

call_one() {
  local sample="$1"
  local depth="$2"
  local bam="$BAM_DIR/${sample}_${depth}.markdup.bam"
  local out_prefix="$VCF_DIR/${sample}_${depth}"

  require_file "$bam"

  gatk --java-options "-Xmx${JAVA_XMX}" Mutect2 \
    -R "$REF_FASTA" \
    -I "$bam" \
    -tumor "$sample" \
    --germline-resource "$GNOMAD_RESOURCE" \
    --panel-of-normals "$PON_RESOURCE" \
    -L "$TARGET_BED" \
    -O "${out_prefix}.unfiltered_TargetRegionOnly.vcf.gz"

  [[ -f "${out_prefix}.unfiltered_TargetRegionOnly.vcf.gz.tbi" ]] || \
    gatk IndexFeatureFile -I "${out_prefix}.unfiltered_TargetRegionOnly.vcf.gz"

  gatk FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "${out_prefix}.unfiltered_TargetRegionOnly.vcf.gz" \
    -O "${out_prefix}.filtered_TargetRegionOnly.vcf.gz"
}

for depth in 5M 20M 100M 200M; do
  call_one mix01 "$depth"
  call_one mix124 "$depth"
done
