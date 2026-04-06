#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# =============================================================================
# 12_group_fragments_by_coordinate.sh
#
# METHODOLOGICAL NOTE — READ BEFORE MODIFYING:
#
# True UMI-based duplicate identification groups reads by the molecular barcode
# ligated to each DNA fragment BEFORE PCR amplification. Because identical
# barcodes prove common origin, UMIs can resolve the rare case of two distinct
# molecules with identical mapping coordinates.
#
# This dataset has NO UMIs. They were stripped during SRA submission:
#   - ENA filereport shows empty submitted_ftp (original files not deposited)
#   - vdb-dump shows FMT: FASTQ (no SPOT_DESC for UMI slots)
#
# Therefore this script implements COORDINATE-BASED CONSENSUS (not UMI calling):
#
#   Groups read pairs by: chr + leftmost_5prime_pos + rightmost_5prime_pos + strand
#
#   This is the SAME grouping assumption as Picard MarkDuplicates. The difference
#   from MarkDuplicates is what happens next: rather than discarding all but the
#   highest-quality duplicate, fgbio calls a CONSENSUS BASE across the family,
#   converting per-read sequencing errors into quality-weighted consensus calls.
#
#   Synthetic 8bp UMIs are derived deterministically from coordinates via MD5:
#     MD5("chr:frag_start:frag_end:strand")[:8] hex nibbles -> ACGT
#     [0-3]->A  [4-7]->C  [8-b]->G  [c-f]->T
#     These UMIs encode NO information beyond the coordinates. They exist solely
#     because fgbio GroupReadsByUmi requires UMI-tagged reads as input.
#     fgbio is therefore invoked with --strategy identity (not adjacency) because
#     edit-distance correction on coordinate-derived UMIs is meaningless.
#
# Limitation vs true UMI calling:
#   Cannot distinguish two distinct molecules with identical 5-prime endpoints.
#   This is rare but non-zero, and represents a systematic undercount of unique
#   molecules. Quantifying the performance gap between this pipeline and a true
#   UMI-aware pipeline motivates preserving UMIs during archival.
#
# Output:
#   ${run}_ontarget_R1/R2.fastq.gz  -- standard FASTQs for Arms A/B (fair baseline)
#   ${run}_umi_R1/R2.fastq.gz       -- UMI-prefixed FASTQs for Arms C/D (fgbio)
#   consensus_metrics/*_family_sizes.tsv
#   consensus_metrics/*_fragment_summary.tsv
# =============================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd bwa
require_cmd samtools
require_cmd python3
require_file "$REF_FASTA"
require_file "$TARGET_BED"

mkdir -p "$RESULTS_DIR/consensus_metrics"

# Write the Python UMI-assignment script to a temp file.
# Reads name-sorted SAM from stdin (piped via samtools view -h).
# Outputs standard and UMI-prefixed FASTQ pairs plus QC metrics.
write_umi_script() {
  local dest="$1"
  cat > "$dest" << 'PYEOF'
"""
Assign deterministic synthetic UMIs to read pairs based on fragment coordinates.

Input:  name-sorted SAM on stdin (samtools view -h pipe)
Output: standard FASTQs (Arms A/B) and UMI-prefixed FASTQs (Arms C/D)

UMI derivation:
  key = "chr:frag_start:frag_end:strand"   (strand = orientation of read 1)
  md5 = MD5(key).hexdigest()
  umi = map each hex nibble [0-3]->A [4-7]->C [8-b]->G [c-f]->T, take first 8
  Deterministic across runs; encodes nothing beyond the coordinates themselves.

fgbio read structure: 8M+T (R1), +T (R2)
  The 8bp UMI is prepended to R1 sequence. R2 is unchanged.
  Synthetic UMI base quality is set to Phred 40 (ASCII 'I').
"""
import sys, gzip, hashlib, collections, argparse

# Map each of the 16 hex characters to one of 4 ACGT bases
HEX_TO_BASE = {}
for _i, _b in enumerate("ACGT"):
    for _h in "0123456789abcdef"[_i * 4 : (_i + 1) * 4]:
        HEX_TO_BASE[_h] = _b

def coord_umi(chrom, start, end, strand):
    key = f"{chrom}:{start}:{end}:{strand}"
    return "".join(HEX_TO_BASE[c] for c in hashlib.md5(key.encode()).hexdigest()[:8])

def write_record(fh, name, seq, qual):
    fh.write(f"@{name}\n{seq}\n+\n{qual}\n".encode())

ap = argparse.ArgumentParser()
ap.add_argument("--out-std-r1", required=True)
ap.add_argument("--out-std-r2", required=True)
ap.add_argument("--out-umi-r1", required=True)
ap.add_argument("--out-umi-r2", required=True)
ap.add_argument("--metrics",    required=True)  # family_size -> group_count TSV
ap.add_argument("--summary",    required=True)  # summary metrics TSV
args = ap.parse_args()

coord_groups = collections.Counter()  # coord_key -> pair count (= family size)
total_pairs  = 0
pending      = {}  # read_name -> (is_read1, is_rev, chrom, pos, seq, qual)
# Note on read names: BWA-MEM does not append /1 or /2 to SAM query names;
# the R1/R2 distinction is encoded in the FLAG field (bits 64/128). The
# pending dict therefore keys on the bare read name from field 0.

with gzip.open(args.out_std_r1, "wb") as s1, \
     gzip.open(args.out_std_r2, "wb") as s2, \
     gzip.open(args.out_umi_r1, "wb") as u1, \
     gzip.open(args.out_umi_r2, "wb") as u2:

    for line in sys.stdin:
        if line.startswith("@"):
            continue  # skip SAM header lines

        f    = line.rstrip("\n").split("\t")
        flag = int(f[1])

        # Skip unmapped, mate-unmapped, secondary, supplementary.
        # -F 2048 in samtools view should have removed supplementary already,
        # but guard here too to prevent pending dict accumulation.
        if flag & 4 or flag & 8 or flag & 256 or flag & 2048:
            continue

        name = f[0]
        seq  = f[9]
        rec  = (bool(flag & 64),   # is_read1
                bool(flag & 16),   # is_reverse
                f[2],              # chrom
                int(f[3]) - 1,     # pos (0-based; SAM is 1-based)
                seq,               # sequence
                f[10])             # base quality string

        if name not in pending:
            pending[name] = rec
        else:
            other = pending.pop(name)
            r1, r2 = (rec, other) if rec[0] else (other, rec)

            # Fragment coordinates.
            # frag_start: leftmost 0-based position across both reads.
            # frag_end:   rightmost end position (pos + read_length) across
            #             both reads. Using pos + len accounts for reads of
            #             different lengths that share a start position —
            #             they represent distinct molecules and must hash
            #             to distinct coordinate groups.
            # strand: orientation of read 1 (same assumption as MarkDuplicates).
            frag_start = min(r1[3], r2[3])
            frag_end   = max(r1[3] + len(r1[4]), r2[3] + len(r2[4]))
            strand     = "-" if r1[1] else "+"

            coord_key = f"{r1[2]}:{frag_start}:{frag_end}:{strand}"
            coord_groups[coord_key] += 1

            umi = coord_umi(r1[2], frag_start, frag_end, strand)
            total_pairs += 1

            seq1, qual1 = r1[4], r1[5]
            seq2, qual2 = r2[4], r2[5]

            # Standard FASTQs -- identical reads, no modification, for Arms A/B.
            # Using the same on-target read set ensures fair comparison across arms.
            write_record(s1, f"{name}/1", seq1, qual1)
            write_record(s2, f"{name}/2", seq2, qual2)

            # UMI-prefixed FASTQs for Arms C/D.
            # fgbio FastqToBam with read structure 8M+T extracts the first 8
            # bases as the UMI into the RX tag. Phred 40 ('I') is used for
            # UMI base quality -- these synthetic bases are high-confidence
            # placeholders that do not affect downstream consensus calling.
            write_record(u1, f"{name}/1", umi + seq1, "I" * 8 + qual1)
            write_record(u2, f"{name}/2", seq2, qual2)

# Safety valve: unpaired reads left in pending were never written.
# This should be zero for a clean name-sorted BAM. Non-zero indicates
# chimeric reads, truncated input, or a name-sort anomaly.
if pending:
    print(f"WARNING: {len(pending)} unpaired reads were not written to output",
          file=sys.stderr)

# Family size distribution: how many coordinate groups have 1, 2, 3... read pairs
family_size_dist = collections.Counter(coord_groups.values())
unique_frags     = len(coord_groups)
dup_rate         = 1.0 - (unique_frags / total_pairs) if total_pairs > 0 else 0.0
mean_family      = total_pairs / unique_frags if unique_frags > 0 else float("nan")

with open(args.metrics, "w") as mf:
    mf.write("family_size\tgroup_count\n")
    for size in sorted(family_size_dist):
        mf.write(f"{size}\t{family_size_dist[size]}\n")

with open(args.summary, "w") as sf:
    sf.write("metric\tvalue\n")
    sf.write(f"total_pairs_processed\t{total_pairs}\n")
    sf.write(f"unique_fragments\t{unique_frags}\n")
    sf.write(f"estimated_duplicate_rate\t{dup_rate:.4f}\n")
    sf.write(f"mean_family_size\t{mean_family:.2f}\n")

print(f"  Total pairs:          {total_pairs:,}", file=sys.stderr)
print(f"  Unique fragments:     {unique_frags:,}", file=sys.stderr)
print(f"  Est. duplicate rate:  {dup_rate:.1%}", file=sys.stderr)
print(f"  Mean family size:     {mean_family:.2f}", file=sys.stderr)
PYEOF
}

process_sample() {
  local run="$1"
  local label="$2"   # "tumor" or "normal"

  local r1="$FASTQ_DIR/${run}_1.fastq.gz"
  local r2="$FASTQ_DIR/${run}_2.fastq.gz"
  require_file "$r1"
  require_file "$r2"

  local out_std_r1="$FASTQ_DIR/${run}_ontarget_R1.fastq.gz"
  local out_std_r2="$FASTQ_DIR/${run}_ontarget_R2.fastq.gz"
  local out_umi_r1="$FASTQ_DIR/${run}_umi_R1.fastq.gz"
  local out_umi_r2="$FASTQ_DIR/${run}_umi_R2.fastq.gz"

  if [[ -f "$out_std_r1" && -f "$out_std_r2" && -f "$out_umi_r1" && -f "$out_umi_r2" ]]; then
    log "[$label] Output FASTQs already exist — skipping."
    return 0
  fi

  local tmp_bam="$TMP_DIR/${run}_tmp_coord.bam"
  local tmp_ontarget="$TMP_DIR/${run}_ontarget_namesort.bam"
  local py_script="$TMP_DIR/assign_umis_${run}_$$.py"

  write_umi_script "$py_script"

  # -------------------------------------------------------------------------
  # Step 1: Temporary coordinate-extraction alignment
  #
  # Aligns the full FASTQ to hg38 to obtain genomic coordinates per read pair.
  # This BAM is temporary — deleted after UMI assignment. Pipeline BAMs for
  # Arms A/B and C/D are produced separately in scripts 13 and 14.
  # -------------------------------------------------------------------------
  if [[ ! -f "${tmp_bam}.bai" ]]; then
    log "[$label] Step 1: Aligning to hg38 for coordinate extraction..."
    bwa mem -t "$THREADS" \
      -R "@RG\tID:${run}\tSM:${label}\tPL:ILLUMINA\tLB:lib1" \
      "$REF_FASTA" "$r1" "$r2" \
      | samtools sort -@ "$THREADS" -T "$TMP_DIR/${run}_sort" -o "$tmp_bam"
    samtools index "$tmp_bam"
    log "[$label] Temporary alignment complete."
  else
    log "[$label] Temporary BAM exists — skipping alignment."
  fi

  # -------------------------------------------------------------------------
  # Step 2: Filter to on-target read pairs overlapping the BRP panel BED
  #
  # -f 1  : read is paired
  # -F 12 : neither read nor mate is unmapped (bit 4 OR bit 8)
  # -L    : read position overlaps target BED
  #
  # NOTE: -L checks only the read's own alignment position, not its mate's.
  # A pair where R1 is on-target and R2 is off-target (or vice versa) will
  # be included based on the on-target read. For amplicon panel sequencing
  # (like BRP), where reads are tightly enriched to target regions, this is
  # acceptable — the vast majority of pairs will have both mates on-target.
  # For hybrid capture, where insert sizes are larger and mate positions more
  # variable, a two-pass approach (extract names, then fetch both mates) would
  # be more correct. Off-target mates will appear as orphans in the name-sorted
  # BAM and are silently skipped by the Python pending-dict logic.
  #
  # Name-sort so the Python script can process each pair together.
  # -------------------------------------------------------------------------
  if [[ ! -f "$tmp_ontarget" ]]; then
    log "[$label] Step 2: Filtering to on-target read pairs..."
    samtools view -b -L "$TARGET_BED" -f 1 -F 12 "$tmp_bam" \
      | samtools sort -n -@ "$THREADS" -T "$TMP_DIR/${run}_nsort" -o "$tmp_ontarget"
    log "[$label] On-target filtering complete."
  else
    log "[$label] On-target BAM exists — skipping filter."
  fi

  # -------------------------------------------------------------------------
  # Steps 3-6: Synthetic UMI assignment and FASTQ output
  #
  # samtools view -h pipes header + alignments to Python, which:
  #   3. Groups pairs by read name (name-sorted BAM guarantees adjacency)
  #   4. Computes fragment coordinates, derives deterministic 8bp UMI via MD5
  #   5. Writes UMI-prefixed FASTQs for Arms C/D (fgbio read structure 8M+T +T)
  #   6. Writes standard FASTQs for Arms A/B (same reads, no modification)
  #      Emits family size distribution and summary QC metrics
  # -------------------------------------------------------------------------
  log "[$label] Steps 3-6: Assigning UMIs and writing FASTQs..."
  samtools view -h "$tmp_ontarget" \
    | python3 "$py_script" \
        --out-std-r1  "$out_std_r1" \
        --out-std-r2  "$out_std_r2" \
        --out-umi-r1  "$out_umi_r1" \
        --out-umi-r2  "$out_umi_r2" \
        --metrics     "$RESULTS_DIR/consensus_metrics/${label}_family_sizes.tsv" \
        --summary     "$RESULTS_DIR/consensus_metrics/${label}_fragment_summary.tsv"

  log "[$label] QC metrics written to $RESULTS_DIR/consensus_metrics/"

  # Temporary alignment BAMs are large — remove them now
  log "[$label] Cleaning up temporary files..."
  rm -f "$tmp_bam" "${tmp_bam}.bai" "$tmp_ontarget" "$py_script"

  log "[$label] Complete."
  log "  Standard FASTQs (Arms A/B): ${run}_ontarget_R1/R2.fastq.gz"
  log "  UMI FASTQs     (Arms C/D): ${run}_umi_R1/R2.fastq.gz"
}

process_sample "$TUMOR_RUN"  "tumor"
process_sample "$NORMAL_RUN" "normal"

log "Script 12 complete. FASTQs ready for scripts 13 (Arms A/B) and 14 (Arms C/D)."
