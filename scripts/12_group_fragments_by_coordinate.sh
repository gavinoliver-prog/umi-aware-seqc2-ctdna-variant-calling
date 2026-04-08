#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# =============================================================================
# 12_group_fragments_by_coordinate.sh  (v2 — parallel pigz compression)
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
#     key = "chr:frag_start:frag_end:strand"
#     umi = first 8 hex nibbles of MD5(key), each mapped: [0-3]->A [4-7]->C [8-b]->G [c-f]->T
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
# PERFORMANCE DESIGN (v2):
#   Original bottleneck: Python's gzip.open() compressing 4 streams serially.
#   Fix: Python writes uncompressed bytes to 4 named pipes (FIFOs); 4 independent
#   pigz processes compress in parallel. Python becomes I/O-bound (fast) rather
#   than CPU-bound on compression (slow). Both samples run in parallel via &.
#
#   Steps 1+2 are now a single pipeline (BWA → samtools view -L → samtools sort -n)
#   eliminating the intermediate coordinate-sorted BAM and index.
#
#   Fallback if still too slow: split namesorted BAM by chromosome prefix and
#   run Python workers in parallel, merging output with cat. Not implemented
#   here because the pigz design should be sufficient on 8+ cores.
#
# PARALLELISM CONTROL:
#   PARALLEL_SAMPLES=1 (default): tumor and normal processed simultaneously.
#     Requires ~14GB RAM free (two concurrent bwa mem processes, ~6-7GB each).
#   PARALLEL_SAMPLES=0: serial execution, lower peak memory.
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
require_cmd pigz
require_cmd python3
require_file "$REF_FASTA"
require_file "$TARGET_BED"

mkdir -p "$RESULTS_DIR/consensus_metrics"

# -----------------------------------------------------------------------------
# RAM availability check
#
# Parallel BWA-MEM (two concurrent processes) requires ~6-7GB each plus
# sort buffers. If available RAM is below 14GB, warn the user and suggest
# running serially via PARALLEL_SAMPLES=0.
# This is advisory only — the script does not abort. The machine may still
# have enough RAM after the OS reclaims cached pages.
# -----------------------------------------------------------------------------
available_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
available_gb=$(( available_kb / 1024 / 1024 ))
if [[ "${PARALLEL_SAMPLES:-1}" == "1" && $available_gb -lt 14 ]]; then
  log "WARNING: Only ~${available_gb}GB RAM currently available."
  log "WARNING: Parallel BWA-MEM (2 samples × ~6-7GB) may exhaust memory and OOM."
  log "WARNING: If the run fails, rerun with PARALLEL_SAMPLES=0 for serial execution."
fi

# Write the Python UMI-assignment script to a temp file.
# Reads name-sorted SAM from stdin (piped via samtools view -h).
# Writes uncompressed FASTQ to the 4 file paths given as arguments.
# Compression is handled externally by pigz reading from named pipes.
write_umi_script() {
  local dest="$1"
  cat > "$dest" << 'PYEOF'
"""
Assign deterministic synthetic UMIs to read pairs based on fragment coordinates.

Input:  name-sorted SAM on stdin (samtools view -h pipe)
Output: uncompressed FASTQ written to 4 file paths (may be named pipes/FIFOs).
        Compression is handled by the caller (pigz).

UMI derivation:
  key = "chr:frag_start:frag_end:strand"   (strand = orientation of read 1)
  md5 = MD5(key).hexdigest()
  umi = map each hex nibble [0-3]->A [4-7]->C [8-b]->G [c-f]->T, take first 8
  Deterministic across runs; encodes nothing beyond the coordinates themselves.

fgbio read structure: 8M+T (R1), +T (R2)
  The 8bp UMI is prepended to R1 sequence. R2 is unchanged.
  Synthetic UMI base quality is set to Phred 40 (ASCII 'I').
"""
import sys, hashlib, collections, argparse

# Map each of the 16 hex characters to one of 4 ACGT bases
HEX_TO_BASE = {}
for _i, _b in enumerate("ACGT"):
    for _h in "0123456789abcdef"[_i * 4 : (_i + 1) * 4]:
        HEX_TO_BASE[_h] = _b

def coord_umi(chrom, start, end, strand):
    key = f"{chrom}:{start}:{end}:{strand}"
    return "".join(HEX_TO_BASE[c] for c in hashlib.md5(key.encode()).hexdigest()[:8])

def write_record(fh, name, seq, qual):
    # fh is a binary file object (open in "wb" mode) pointing to a named pipe.
    # Writing bytes avoids a text-mode encode round-trip and is compatible
    # with pigz reading raw bytes from the other end of the pipe.
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
# the R1/R2 distinction is encoded in the FLAG field (bits 64/128).

# Open all 4 outputs as binary — works for both regular files and named pipes.
# Named pipes: pigz readers are already running before Python opens these;
# open() blocks until the other end (pigz) also opens the pipe, which
# guarantees no data is lost even before Python's first write.
with open(args.out_std_r1, "wb") as s1, \
     open(args.out_std_r2, "wb") as s2, \
     open(args.out_umi_r1, "wb") as u1, \
     open(args.out_umi_r2, "wb") as u2:

    for line in sys.stdin:
        if line.startswith("@"):
            continue  # skip SAM header lines

        f    = line.rstrip("\n").split("\t")
        flag = int(f[1])

        # Skip unmapped, mate-unmapped, secondary, supplementary.
        # samtools view -F 2316 upstream removes most of these, but
        # guard here to prevent pending dict accumulation on edge cases.
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
# Should be zero for a clean name-sorted BAM. Non-zero indicates chimeric
# reads, truncated input, or a name-sort anomaly.
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

  local tmp_ontarget="$TMP_DIR/${run}_ontarget_namesort.bam"
  local py_script="$TMP_DIR/assign_umis_${run}_$$.py"
  local pipe_dir="$TMP_DIR/${run}_pipes_$$"

  # Clean up named pipes and partial outputs on any failure in this function
  cleanup() {
    local exit_code=$?
    [[ $exit_code -ne 0 ]] && log "[$label] ERROR: cleaning up partial outputs (exit $exit_code)"
    rm -f "${out_std_r1}.tmp" "${out_std_r2}.tmp" "${out_umi_r1}.tmp" "${out_umi_r2}.tmp"
    rm -rf "$pipe_dir"
  }
  trap cleanup EXIT

  write_umi_script "$py_script"

  # -------------------------------------------------------------------------
  # Steps 1+2 combined: Align → filter on-target → name-sort
  #
  # v1 used two separate steps: coord-sort+index (for samtools view -L),
  # then view -L | sort -n. This required writing and reading a full
  # coordinate-sorted BAM (~30-40GB for these samples).
  #
  # v2 streams directly: BWA → samtools view -L → samtools sort -n.
  # samtools view -L on a stream scans all records without an index.
  # For panel sequencing (BRP: ~90% on-target rate), the filter eliminates
  # most reads before they hit the sort buffer, reducing sort memory and I/O.
  #
  # -F 2316: skip unmapped(4) + mate-unmapped(8) + secondary(256) + suppl(2048)
  # -------------------------------------------------------------------------
  if [[ ! -f "$tmp_ontarget" ]]; then
    log "[$label] Steps 1-2: Align → filter on-target → name-sort (single pipeline)..."
    bwa mem -t "$THREADS" \
      -R "@RG\tID:${run}\tSM:${label}\tPL:ILLUMINA\tLB:lib1" \
      "$REF_FASTA" "$r1" "$r2" \
      | samtools view -u -f 1 -F 2316 -L "$TARGET_BED" \
      | samtools sort -n -@ "$THREADS" -m 4G \
                      -T "$TMP_DIR/${run}_nsort" \
                      -o "$tmp_ontarget"
    log "[$label] Name-sorted on-target BAM ready."
  else
    log "[$label] Name-sorted BAM exists — skipping alignment."
  fi

  # -------------------------------------------------------------------------
  # Steps 3-6: UMI assignment with named-pipe + pigz parallel compression
  #
  # Architecture:
  #   [samtools view] → [Python] → [FIFO s1] → [pigz] → out_std_r1.fastq.gz
  #                              → [FIFO s2] → [pigz] → out_std_r2.fastq.gz
  #                              → [FIFO u1] → [pigz] → out_umi_r1.fastq.gz
  #                              → [FIFO u2] → [pigz] → out_umi_r2.fastq.gz
  #
  # Python writes uncompressed bytes to each FIFO. Each pigz process reads
  # from its own FIFO and compresses independently. The 4 pigz processes run
  # concurrently so compression throughput is 4× that of a single gzip.
  #
  # pigz thread count: THREADS/4, floored to 1, capped at 4. When two
  # process_sample calls run in parallel (tumor + normal), the combined
  # pigz load is 2×(4 streams × pigz_t). Capping at 4 prevents the 8-core
  # machine from being saturated by compression alone.
  # -------------------------------------------------------------------------
  log "[$label] Steps 3-6: UMI assignment + parallel pigz compression..."

  mkdir -p "$pipe_dir"
  local p_std_r1="$pipe_dir/std_r1"
  local p_std_r2="$pipe_dir/std_r2"
  local p_umi_r1="$pipe_dir/umi_r1"
  local p_umi_r2="$pipe_dir/umi_r2"
  mkfifo "$p_std_r1" "$p_std_r2" "$p_umi_r1" "$p_umi_r2"

  local pigz_t=$(( THREADS / 4 ))
  [[ $pigz_t -lt 1 ]] && pigz_t=1
  [[ $pigz_t -gt 4 ]] && pigz_t=4

  # Start 4 pigz compressors in background, each dedicated to one FIFO.
  # Write to .tmp paths; renamed to final names only on full success to
  # prevent a partial compressed file from passing the idempotency check.
  pigz -p "$pigz_t" < "$p_std_r1" > "${out_std_r1}.tmp" & local pid_s1=$!
  pigz -p "$pigz_t" < "$p_std_r2" > "${out_std_r2}.tmp" & local pid_s2=$!
  pigz -p "$pigz_t" < "$p_umi_r1" > "${out_umi_r1}.tmp" & local pid_u1=$!
  pigz -p "$pigz_t" < "$p_umi_r2" > "${out_umi_r2}.tmp" & local pid_u2=$!

  # Python opens each FIFO (blocks until the pigz reader above has opened
  # the other end — guaranteed since all 4 pigz processes are already running).
  # On exit, Python closes all 4 FIFOs; each pigz reader reaches EOF, flushes
  # its compressed output, and exits cleanly.
  samtools view -h "$tmp_ontarget" \
    | python3 "$py_script" \
        --out-std-r1  "$p_std_r1" \
        --out-std-r2  "$p_std_r2" \
        --out-umi-r1  "$p_umi_r1" \
        --out-umi-r2  "$p_umi_r2" \
        --metrics     "$RESULTS_DIR/consensus_metrics/${label}_family_sizes.tsv" \
        --summary     "$RESULTS_DIR/consensus_metrics/${label}_fragment_summary.tsv"

  # Wait for all pigz processes; any non-zero exit means corrupted output
  local fail=0
  wait "$pid_s1" || { log "[$label] ERROR: pigz failed on std_r1"; fail=1; }
  wait "$pid_s2" || { log "[$label] ERROR: pigz failed on std_r2"; fail=1; }
  wait "$pid_u1" || { log "[$label] ERROR: pigz failed on umi_r1"; fail=1; }
  wait "$pid_u2" || { log "[$label] ERROR: pigz failed on umi_r2"; fail=1; }

  if [[ $fail -ne 0 ]]; then
    log "[$label] Compression failed — partial outputs will be removed by cleanup trap"
    return 1
  fi

  # Atomically rename .tmp files to final names only after all 4 succeed
  mv "${out_std_r1}.tmp" "$out_std_r1"
  mv "${out_std_r2}.tmp" "$out_std_r2"
  mv "${out_umi_r1}.tmp" "$out_umi_r1"
  mv "${out_umi_r2}.tmp" "$out_umi_r2"

  rm -rf "$pipe_dir"
  log "[$label] QC metrics written to $RESULTS_DIR/consensus_metrics/"

  log "[$label] Cleaning up temporary alignment BAM..."
  rm -f "$tmp_ontarget" "$py_script"

  # Disarm the cleanup trap — outputs are complete
  trap - EXIT

  log "[$label] Complete."
  log "  Standard FASTQs (Arms A/B): ${run}_ontarget_R1/R2.fastq.gz"
  log "  UMI FASTQs     (Arms C/D): ${run}_umi_R1/R2.fastq.gz"
}

# =============================================================================
# Run tumor and normal samples — parallel (default) or serial.
#
# PARALLEL_SAMPLES=1 (default):
#   Both samples processed simultaneously. Peak CPU = 2×THREADS during
#   alignment; peak RAM = ~14GB (two concurrent bwa mem + sort buffers).
#   Total wall time is lower because alignment and Python/pigz phases of
#   each sample overlap with the other sample's processing.
#
# PARALLEL_SAMPLES=0:
#   Serial execution. Suitable for memory-constrained systems (<14GB free)
#   or for cleaner log output. Set in environment or prefix the command:
#     PARALLEL_SAMPLES=0 bash scripts/12_group_fragments_by_coordinate.sh
# =============================================================================
log "PARALLEL_SAMPLES=${PARALLEL_SAMPLES:-1} (available RAM: ~${available_gb}GB)"

fail=0
if [[ "${PARALLEL_SAMPLES:-1}" == "1" ]]; then
  process_sample "$TUMOR_RUN"  "tumor"  & tumor_pid=$!
  process_sample "$NORMAL_RUN" "normal" & normal_pid=$!
  wait "$tumor_pid"  || { log "ERROR: tumor sample processing failed";  fail=1; }
  wait "$normal_pid" || { log "ERROR: normal sample processing failed"; fail=1; }
else
  log "Running serially (PARALLEL_SAMPLES=0)..."
  process_sample "$TUMOR_RUN"  "tumor"  || fail=1
  process_sample "$NORMAL_RUN" "normal" || fail=1
fi

if [[ $fail -ne 0 ]]; then
  log "ERROR: one or more samples failed — check logs above"
  exit 1
fi

log "Script 12 complete. FASTQs ready for scripts 13 (Arms A/B) and 14 (Arms C/D)."
