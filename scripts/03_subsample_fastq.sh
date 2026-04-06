#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

require_cmd seqtk
require_cmd zcat
require_cmd gzip

cd "$FASTQ_DIR"

for run in "$MIX01_RUN" "$MIX124_RUN"; do
  seqtk sample -s100 "${run}_1.fastq.gz" 5000000  | gzip > "sub_5M_${run}_R1.fastq.gz"
  seqtk sample -s100 "${run}_2.fastq.gz" 5000000  | gzip > "sub_5M_${run}_R2.fastq.gz"

  seqtk sample -s100 "${run}_1.fastq.gz" 20000000 | gzip > "sub_20M_${run}_R1.fastq.gz"
  seqtk sample -s100 "${run}_2.fastq.gz" 20000000 | gzip > "sub_20M_${run}_R2.fastq.gz"
done

# Large random subsets were impractical with seqtk due to memory use at this scale.
for run in "$MIX01_RUN" "$MIX124_RUN"; do
  zcat "${run}_1.fastq.gz" | head -400000000 | gzip > "sub_100M_${run}_R1.fastq.gz"
  zcat "${run}_2.fastq.gz" | head -400000000 | gzip > "sub_100M_${run}_R2.fastq.gz"

  zcat "${run}_1.fastq.gz" | head -800000000 | gzip > "sub_200M_${run}_R1.fastq.gz"
  zcat "${run}_2.fastq.gz" | head -800000000 | gzip > "sub_200M_${run}_R2.fastq.gz"
done
