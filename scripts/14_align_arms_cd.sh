#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Source the fgbio JAR path written by 01_setup_env.sh
FGBIO_ENV="$(dirname "$0")/fgbio_env.sh"
[[ -f "$FGBIO_ENV" ]] && source "$FGBIO_ENV"
[[ -n "${FGBIO_JAR:-}" ]] || { echo "ERROR: FGBIO_JAR not set. Run 01_setup_env.sh first." >&2; exit 1; }
require_file "$FGBIO_JAR"

# =============================================================================
# Consensus pipeline for Arms C and D.
#
# Input:  UMI-prefixed FASTQs from script 12 (${run}_umi_R1/R2.fastq.gz)
#         Read structure: 8M+T (R1) — 8bp synthetic UMI + template
#                         +T   (R2) — template only
#
# Pipeline:
#   1. fgbio FastqToBam    — extract 8bp UMI from R1 prefix into RX tag
#   2. BWA-MEM             — align reads; preserve RX tag via fgbio ZipperBams
#   3. fgbio GroupReadsByUmi --strategy identity
#                          — group reads with identical UMIs into families
#                            identity (not adjacency): UMIs are coordinate-
#                            derived, so edit-distance correction is meaningless
#   4. fgbio CallMolecularConsensusReads --min-reads $UMI_MIN_READS
#                          — collapse each family into a consensus read pair;
#                            base quality reflects agreement across family members
#   5. fgbio FilterConsensusReads --min-base-quality $UMI_MIN_BASE_QUALITY
#                          — mask low-confidence consensus bases and remove
#                            families that fall below quality thresholds
#   6. BWA-MEM re-alignment — consensus reads are unmapped after step 5;
#                             re-align to get final coordinate-sorted BAM
#
# Arms C and D share the same consensus tumor BAM (same as A/B share markdup BAM).
# Output: ARM_C_DIR/tumor.consensus.bam
#         ARM_D_DIR/tumor.consensus.bam  (symlink to ARM_C_DIR version)
#         ARM_D_DIR/normal.consensus.bam
# =============================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd bwa
require_cmd samtools
require_cmd gatk
require_file "$REF_FASTA"

mkdir -p "$ARM_C_DIR" "$ARM_D_DIR"

process_sample() {
  local run="$1"
  local label="$2"    # "tumor" or "normal"; must match Mutect2 -tumor/-normal in script 16
  local out_bam="$3"

  local umi_r1="$FASTQ_DIR/${run}_umi_R1.fastq.gz"
  local umi_r2="$FASTQ_DIR/${run}_umi_R2.fastq.gz"
  require_file "$umi_r1"
  require_file "$umi_r2"

  local out_dir
  out_dir="$(dirname "$out_bam")"

  if [[ -f "${out_bam}.bai" ]]; then
    log "[$label] $out_bam already indexed — skipping."
    return 0
  fi

  # Intermediate BAMs in TMP_DIR — all cleaned up at end
  local ubam="$TMP_DIR/${run}_cd_ubam.bam"
  local mapped_nsort="$TMP_DIR/${run}_cd_mapped_nsort.bam"
  local merged_bam="$TMP_DIR/${run}_cd_merged.bam"
  local grouped_bam="$TMP_DIR/${run}_cd_grouped.bam"
  local consensus_raw="$TMP_DIR/${run}_cd_consensus_raw.bam"
  local consensus_filt="$TMP_DIR/${run}_cd_consensus_filt.bam"
  local consensus_r1="$TMP_DIR/${run}_cd_consensus_R1.fastq.gz"
  local consensus_r2="$TMP_DIR/${run}_cd_consensus_R2.fastq.gz"

  # -------------------------------------------------------------------------
  # Step 1: FastqToBam — extract 8bp synthetic UMI from R1 prefix into RX tag
  #
  # Read structure "8M+T" means: 8 molecular barcode bases (M) then remaining
  # template bases (T). "+T" for R2 means all bases are template.
  # fgbio moves the M bases into the RX BAM tag and trims them from the read
  # sequence, leaving clean template reads for alignment.
  # -------------------------------------------------------------------------
  if [[ ! -f "$ubam" ]]; then
    log "[$label] Step 1: FastqToBam — extracting UMI tags..."
    java -Xmx"${JAVA_XMX}" -jar "$FGBIO_JAR" FastqToBam \
      --input "$umi_r1" "$umi_r2" \
      --output "$ubam" \
      --read-structures "8M+T" "+T" \
      --sample "$label" \
      --library "lib1" \
      --read-group-id "$run" \
      --platform ILLUMINA \
      --sort true
  else
    log "[$label] Step 1: uBAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 2a: BWA-MEM alignment
  #
  # samtools fastq converts the uBAM back to FASTQ for BWA-MEM.
  # The aligned output is name-sorted for ZipperBams in step 2b.
  # We do NOT use -C (pass FASTQ comments) here — ZipperBams handles tag
  # transfer from the uBAM, so the mapped BAM only needs alignment coordinates.
  # -------------------------------------------------------------------------
  if [[ ! -f "$mapped_nsort" ]]; then
    log "[$label] Step 2a: BWA-MEM alignment..."
    samtools fastq -N "$ubam" \
      | bwa mem -t "$THREADS" -K 100000000 -p \
          -R "@RG\tID:${run}\tSM:${label}\tPL:ILLUMINA\tLB:lib1" \
          "$REF_FASTA" - \
      | samtools sort -n -@ "$THREADS" -T "$TMP_DIR/${run}_cd_nsort" -o "$mapped_nsort"
  else
    log "[$label] Step 2a: Mapped BAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 2b: ZipperBams — merge UMI tags from uBAM into aligned BAM
  #
  # BWA-MEM discards BAM tags from input. ZipperBams re-attaches tags from
  # the original uBAM (including RX) to the aligned records by matching
  # read names. Both inputs must be in the same query-name order.
  # Output is query-name sorted; we coordinate-sort afterward for GRBU.
  # -------------------------------------------------------------------------
  if [[ ! -f "$merged_bam" ]]; then
    log "[$label] Step 2b: ZipperBams — merging UMI tags into aligned BAM..."
    java -Xmx"${JAVA_XMX}" -jar "$FGBIO_JAR" ZipperBams \
      --unmapped "$ubam" \
      --mapped "$mapped_nsort" \
      --output "$merged_bam" \
      --tags RX
    # Coordinate-sort for GroupReadsByUmi
    samtools sort -@ "$THREADS" -T "$TMP_DIR/${run}_cd_coord" -o "${merged_bam%.bam}.coord.bam" "$merged_bam"
    mv "${merged_bam%.bam}.coord.bam" "$merged_bam"
    samtools index "$merged_bam"
    rm -f "$mapped_nsort"
  else
    log "[$label] Step 2b: Merged BAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 3: GroupReadsByUmi --strategy identity
  #
  # Groups read pairs that share an identical RX tag (synthetic UMI) into
  # a molecule family, tagged with MI (molecule identifier).
  #
  # --strategy identity: exact UMI match only. Adjacency/edit-distance
  # correction is NOT used because our UMIs are coordinate-derived — two
  # reads with a 1-base UMI difference actually came from different coordinate
  # groups and represent distinct molecules, not PCR errors in the barcode.
  # -------------------------------------------------------------------------
  if [[ ! -f "$grouped_bam" ]]; then
    log "[$label] Step 3: GroupReadsByUmi (strategy=identity)..."
    java -Xmx"${JAVA_XMX}" -jar "$FGBIO_JAR" GroupReadsByUmi \
      --input "$merged_bam" \
      --output "$grouped_bam" \
      --strategy identity \
      --family-size-histogram "${out_dir}/${label}.grbu_family_sizes.txt"
    rm -f "$merged_bam" "${merged_bam}.bai"
  else
    log "[$label] Step 3: Grouped BAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 4: CallMolecularConsensusReads
  #
  # Collapses each read family into a single consensus read pair. For each
  # base position, the consensus base and quality are derived from a
  # quality-weighted vote across all family members. Bases that disagree
  # across family members receive reduced consensus quality — this is the
  # core error-suppression mechanism that makes consensus calling superior
  # to simple duplicate discarding for low-VAF variant detection.
  #
  # --min-reads: families with fewer than this many reads are discarded.
  #   Value 2 removes singletons, which are most likely PCR or sequencing
  #   artefacts rather than true molecules.
  # -------------------------------------------------------------------------
  if [[ ! -f "$consensus_raw" ]]; then
    log "[$label] Step 4: CallMolecularConsensusReads (min-reads=$UMI_MIN_READS)..."
    java -Xmx"${JAVA_XMX}" -jar "$FGBIO_JAR" CallMolecularConsensusReads \
      --input "$grouped_bam" \
      --output "$consensus_raw" \
      --min-reads "$UMI_MIN_READS" \
      --read-group-id "$run"
    rm -f "$grouped_bam"
  else
    log "[$label] Step 4: Raw consensus BAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 5: FilterConsensusReads
  #
  # Applies quality filters to the consensus BAM:
  # --min-base-quality: mask consensus bases below this Phred score to N.
  #   Bases masked to N are not used in variant calling.
  # --min-reads: redundant with step 4 but ensures consistency.
  #
  # Output is an unmapped BAM — consensus reads have no alignment coordinates.
  # -------------------------------------------------------------------------
  if [[ ! -f "$consensus_filt" ]]; then
    log "[$label] Step 5: FilterConsensusReads (min-base-quality=$UMI_MIN_BASE_QUALITY)..."
    java -Xmx"${JAVA_XMX}" -jar "$FGBIO_JAR" FilterConsensusReads \
      --input "$consensus_raw" \
      --output "$consensus_filt" \
      --ref "$REF_FASTA" \
      --min-reads "$UMI_MIN_READS" \
      --min-base-quality "$UMI_MIN_BASE_QUALITY"
    rm -f "$consensus_raw"
  else
    log "[$label] Step 5: Filtered consensus BAM exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # Step 6: Re-align consensus reads
  #
  # Consensus reads from FilterConsensusReads are unmapped. Re-align to hg38
  # to produce the final coordinate-sorted BAM for variant calling.
  # Consensus reads are shorter than input reads (family size >= 2 required)
  # and have improved base qualities — Mutect2 and LoFreq will use these
  # qualities directly in their statistical models.
  #
  # -N: add /1 /2 suffixes to read names for BWA paired-end handling
  # -F 0x900: exclude secondary (0x100) and supplementary (0x800) alignments
  #   — consensus BAMs should not contain these, but guard defensively
  # -------------------------------------------------------------------------
  log "[$label] Step 6: Re-aligning consensus reads..."
  samtools fastq -N -F 0x900 \
    -1 "$consensus_r1" \
    -2 "$consensus_r2" \
    -0 /dev/null -s /dev/null \
    "$consensus_filt"

  bwa mem -t "$THREADS" -K 100000000 \
    -R "@RG\tID:${run}\tSM:${label}\tPL:ILLUMINA\tLB:lib1" \
    "$REF_FASTA" "$consensus_r1" "$consensus_r2" \
    | samtools sort -@ "$THREADS" -T "$TMP_DIR/${run}_cd_final_sort" -o "$out_bam"

  samtools index "$out_bam"

  # Collect post-consensus alignment metrics
  log "[$label] Collecting metrics..."
  samtools flagstat -@ "$THREADS" "$out_bam" \
    > "${out_dir}/${label}.flagstat.txt"

  # CollectInsertSizeMetrics can fail on low read count BAMs (e.g. if consensus
  # calling was very aggressive) — treat as non-fatal
  gatk --java-options "-Xmx${JAVA_XMX}" CollectInsertSizeMetrics \
    -I "$out_bam" \
    -O "${out_dir}/${label}.insert_size_metrics.txt" \
    -H "${out_dir}/${label}.insert_size_histogram.pdf" \
    || log "[$label] CollectInsertSizeMetrics failed — skipping (insufficient reads for histogram)"

  # Clean up remaining temp files
  log "[$label] Cleaning up temporary files..."
  rm -f "$ubam" "$consensus_filt" "$consensus_r1" "$consensus_r2"

  log "[$label] Done: $out_bam"
}

process_sample "$TUMOR_RUN"  "tumor"  "$ARM_C_DIR/tumor.consensus.bam"

# Arms C and D share the same consensus tumor BAM — only calling differs.
# Create a symlink in ARM_D_DIR so script 16 can reference it uniformly.
if [[ ! -L "$ARM_D_DIR/tumor.consensus.bam" ]]; then
  ln -sf "$ARM_C_DIR/tumor.consensus.bam"     "$ARM_D_DIR/tumor.consensus.bam"
  ln -sf "$ARM_C_DIR/tumor.consensus.bam.bai" "$ARM_D_DIR/tumor.consensus.bam.bai"
  log "Symlinked tumor consensus BAM into $ARM_D_DIR"
fi

process_sample "$NORMAL_RUN" "normal" "$ARM_D_DIR/normal.consensus.bam"

log "Script 14 complete."
log "  Consensus tumor BAM  (Arms C+D): $ARM_C_DIR/tumor.consensus.bam"
log "  Consensus normal BAM (Arm D):    $ARM_D_DIR/normal.consensus.bam"
