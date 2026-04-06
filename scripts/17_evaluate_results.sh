#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Core evaluation against the SEQC2 truth set.
# Evaluates all four Mutect2 arms and LoFreq for Arms A and C.
#
# Outputs:
#   evaluation/per_variant_calls.tsv    — one row per known positive, all arms
#   evaluation/per_arm_summary.tsv      — TP/FP/FN/sensitivity/F1 per arm
#   evaluation/cross_arm_comparison.tsv — four key pairwise arm comparisons
#   evaluation/fp_rate_comparison.tsv   — FP rate per Mb per arm

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd bcftools
require_cmd python3
require_file "$KNOWN_POSITIVES_VCF"
require_file "$KNOWN_NEGATIVES_BED"

EVAL_DIR="$RESULTS_DIR/evaluation"
mkdir -p "$EVAL_DIR"

# Known-negative space in Mb (53,185 bp overlapping BRP panel — from SEQC2)
KNOWN_NEG_MB="0.053185"

# Ensure truth VCF is tabix-indexed — required by bcftools isec and view
if [[ ! -f "${KNOWN_POSITIVES_VCF}.tbi" ]]; then
  log "Indexing truth VCF..."
  bcftools index -t "$KNOWN_POSITIVES_VCF"
fi

if [[ -f "$EVAL_DIR/per_variant_calls.tsv"    && \
      -f "$EVAL_DIR/per_arm_summary.tsv"       && \
      -f "$EVAL_DIR/cross_arm_comparison.tsv"  && \
      -f "$EVAL_DIR/fp_rate_comparison.tsv" ]]; then
  log "All evaluation outputs exist — skipping."
  exit 0
fi

# -------------------------------------------------------------------------
# extract_pass_tsv: PASS variants from Mutect2 filtered VCF.
# Uses -s tumor to extract tumor sample's FORMAT/AF, which handles both
# tumor-only (single sample) and tumor-normal (two sample) VCFs uniformly.
# -------------------------------------------------------------------------
extract_pass_tsv() {
  local vcf="$1"
  local out="$2"
  [[ -f "$out" ]] && return 0
  require_file "$vcf"
  bcftools view -f PASS "$vcf" \
    | bcftools query -s tumor -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' \
    > "$out"
}

# -------------------------------------------------------------------------
# extract_lofreq_tsv: variants from LoFreq VCF.
# LoFreq stores AF in INFO (not FORMAT). Uses -f 'PASS,.' to include both
# explicitly-PASS variants and variants with FILTER='.' (produced when
# --no-default-filter is used in script 16).
# -------------------------------------------------------------------------
extract_lofreq_tsv() {
  local vcf="$1"
  local out="$2"
  [[ -f "$out" ]] && return 0
  if [[ ! -f "$vcf" ]]; then
    log "  WARNING: $vcf not found — writing empty placeholder."
    > "$out"
    return 0
  fi
  bcftools view -f 'PASS,.' "$vcf" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' \
    > "$out"
}

# -------------------------------------------------------------------------
# count_fps: PASS variants overlapping known-negative regions.
# Uses -T (targets-file) to avoid tabix requirement on .bed.gz.
# Uses an if/else block rather than a variable holding a flag string, to
# avoid quoting issues when bcftools expands 'PASS,.' from a variable.
# -------------------------------------------------------------------------
count_fps() {
  local vcf="$1"
  local is_lofreq="${2:-0}"
  [[ -f "$vcf" ]] || { echo 0; return 0; }
  if [[ "$is_lofreq" == "1" ]]; then
    bcftools view -f 'PASS,.' -T "$KNOWN_NEGATIVES_BED" "$vcf" \
      | bcftools view -H | wc -l
  else
    bcftools view -f PASS -T "$KNOWN_NEGATIVES_BED" "$vcf" \
      | bcftools view -H | wc -l
  fi
}

# -------------------------------------------------------------------------
# Step 1: Extract PASS call TSVs for all arms and callers
# -------------------------------------------------------------------------
log "Step 1: Extracting PASS variant TSVs..."
extract_pass_tsv "$ARM_A_DIR/mutect2.filtered.vcf.gz" "$EVAL_DIR/arm_a_pass.tsv"
extract_pass_tsv "$ARM_B_DIR/mutect2.filtered.vcf.gz" "$EVAL_DIR/arm_b_pass.tsv"
extract_pass_tsv "$ARM_C_DIR/mutect2.filtered.vcf.gz" "$EVAL_DIR/arm_c_pass.tsv"
extract_pass_tsv "$ARM_D_DIR/mutect2.filtered.vcf.gz" "$EVAL_DIR/arm_d_pass.tsv"
extract_lofreq_tsv "$ARM_A_DIR/lofreq.vcf.gz"         "$EVAL_DIR/arm_a_lofreq_pass.tsv"
extract_lofreq_tsv "$ARM_C_DIR/lofreq.vcf.gz"         "$EVAL_DIR/arm_c_lofreq_pass.tsv"

# -------------------------------------------------------------------------
# Step 2: Count FPs in known-negative space per arm
# -------------------------------------------------------------------------
log "Step 2: Counting false positives in known-negative space..."
fp_arm_a=$(count_fps        "$ARM_A_DIR/mutect2.filtered.vcf.gz")
fp_arm_b=$(count_fps        "$ARM_B_DIR/mutect2.filtered.vcf.gz")
fp_arm_c=$(count_fps        "$ARM_C_DIR/mutect2.filtered.vcf.gz")
fp_arm_d=$(count_fps        "$ARM_D_DIR/mutect2.filtered.vcf.gz")
fp_arm_a_lofreq=$(count_fps "$ARM_A_DIR/lofreq.vcf.gz" 1)
fp_arm_c_lofreq=$(count_fps "$ARM_C_DIR/lofreq.vcf.gz" 1)
log "  FP counts — A:$fp_arm_a B:$fp_arm_b C:$fp_arm_c D:$fp_arm_d A_lofreq:$fp_arm_a_lofreq C_lofreq:$fp_arm_c_lofreq"

# -------------------------------------------------------------------------
# Step 3: Python — build all evaluation TSVs
# -------------------------------------------------------------------------
log "Step 3: Building evaluation TSVs..."

PY_SCRIPT="$TMP_DIR/evaluate_$$.py"
cat > "$PY_SCRIPT" << 'PYEOF'
"""
Evaluation engine for the four-arm ctDNA pipeline benchmark.

Reads truth VCF and per-arm PASS call TSVs (produced by bcftools query).
Produces four output TSVs covering variant-level calls, arm-level summary,
cross-arm comparisons, and FP rate table.
"""
import sys, os, gzip, csv, argparse, collections

DILUTION        = 25.0    # Ef = Sample A diluted 1:24 into Sample B; VAF / 25
NEG_MB_DEFAULT  = 0.053185

def classify_tier(vaf):
    """Tier based on Sample A VAF (pre-dilution)."""
    if vaf >= 0.50: return 1   # >2% expected in Ef
    if vaf >= 0.20: return 2   # 0.8-2%
    if vaf >= 0.05: return 3   # 0.2-0.8%
    return 4                   # <0.2%

TIER_RANGES = {1: '>=50%', 2: '20-50%', 3: '5-20%', 4: '<5%'}

def open_file(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def read_truth_vcf(path):
    """
    Parse known-positives VCF. INFO/VAF = Sample A VAF (per CLAUDE.md).
    Falls back to INFO/AF if VAF not present.
    Tries INFO/GENE, INFO/Gene, INFO/SYMBOL, SnpEff ANN for gene name.

    Note on ALT splitting: split(',')[0] takes the first allele for
    multi-allelic records. MNVs (multi-nucleotide variants) are represented
    as single ALT entries in this truth set (e.g. "AT" not "A,T"), so
    split(',')[0] is a no-op for MNVs and handles true multi-allelics
    by taking the primary allele.

    Returns OrderedDict preserving VCF order: (chrom, pos, ref, alt) -> dict.
    """
    truth = collections.OrderedDict()
    with open_file(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            chrom, pos, ref = p[0], int(p[1]), p[3]
            alt = p[4].split(',')[0]
            info = {}
            for kv in p[7].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    info[k] = v
            vaf = float(info.get('VAF', info.get('AF', '0')))
            gene = info.get('GENE', info.get('Gene', info.get('SYMBOL', '.')))
            if gene == '.' and 'ANN' in info:
                # SnpEff ANN: allele|effect|impact|gene|...
                ann_parts = info['ANN'].split('|')
                gene = ann_parts[3] if len(ann_parts) > 3 else '.'
            truth[(chrom, pos, ref, alt)] = {'vaf': vaf, 'gene': gene}
    return truth

def read_pass_tsv(path):
    """
    Read bcftools query TSV: CHROM POS REF ALT AF.
    Returns dict: (chrom, pos, ref, alt) -> af_string.
    Handles multi-allelic ALT and AF by taking the first value.
    """
    calls = {}
    if not path or not os.path.exists(path):
        return calls
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            p = line.split('\t')
            if len(p) < 5:
                continue
            chrom, pos, ref = p[0], int(p[1]), p[2]
            alt = p[3].split(',')[0]
            af  = p[4].split(',')[0] if p[4] not in ('.', '') else '.'
            calls[(chrom, pos, ref, alt)] = af
    return calls

def f1_precision_recall(tp, fp, fn):
    """
    Compute F1, precision, and recall.

    Note: fp here is the count of calls in the known-negative space only
    (53,185 bp). True FPs outside that space are unknown. Therefore precision
    is a LOWER-BOUND ESTIMATE — the true precision may be higher if the
    unknown-space calls are mostly true positives (unlikely at low VAF, but
    possible). Use with appropriate caution in portfolio write-up.
    """
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return f1, prec, rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth-vcf',       required=True)
    ap.add_argument('--arm-a',           required=True)
    ap.add_argument('--arm-b',           required=True)
    ap.add_argument('--arm-c',           required=True)
    ap.add_argument('--arm-d',           required=True)
    ap.add_argument('--arm-a-lofreq',    default='')
    ap.add_argument('--arm-c-lofreq',    default='')
    ap.add_argument('--fp-arm-a',        type=int, default=0)
    ap.add_argument('--fp-arm-b',        type=int, default=0)
    ap.add_argument('--fp-arm-c',        type=int, default=0)
    ap.add_argument('--fp-arm-d',        type=int, default=0)
    ap.add_argument('--fp-arm-a-lofreq', type=int, default=0)
    ap.add_argument('--fp-arm-c-lofreq', type=int, default=0)
    ap.add_argument('--known-neg-mb',    type=float, default=NEG_MB_DEFAULT)
    ap.add_argument('--out-variants',    required=True)
    ap.add_argument('--out-summary',     required=True)
    ap.add_argument('--out-cross-arm',   required=True)
    ap.add_argument('--out-fp-rates',    required=True)
    args = ap.parse_args()

    truth = read_truth_vcf(args.truth_vcf)

    ARMS = ['arm_a', 'arm_b', 'arm_c', 'arm_d', 'arm_a_lofreq', 'arm_c_lofreq']
    calls = {
        'arm_a':        read_pass_tsv(args.arm_a),
        'arm_b':        read_pass_tsv(args.arm_b),
        'arm_c':        read_pass_tsv(args.arm_c),
        'arm_d':        read_pass_tsv(args.arm_d),
        'arm_a_lofreq': read_pass_tsv(args.arm_a_lofreq),
        'arm_c_lofreq': read_pass_tsv(args.arm_c_lofreq),
    }
    fp_counts = {
        'arm_a': args.fp_arm_a, 'arm_b': args.fp_arm_b,
        'arm_c': args.fp_arm_c, 'arm_d': args.fp_arm_d,
        'arm_a_lofreq': args.fp_arm_a_lofreq,
        'arm_c_lofreq': args.fp_arm_c_lofreq,
    }
    arm_labels = {
        'arm_a': ('A', 'Mutect2'), 'arm_b': ('B', 'Mutect2'),
        'arm_c': ('C', 'Mutect2'), 'arm_d': ('D', 'Mutect2'),
        'arm_a_lofreq': ('A', 'LoFreq'), 'arm_c_lofreq': ('C', 'LoFreq'),
    }

    # Accumulate tier-level TP/FN per arm
    tier_stats  = {t: {a: {'tp': 0, 'fn': 0} for a in ARMS} for t in [1,2,3,4]}
    tier_totals = {1: 0, 2: 0, 3: 0, 4: 0}

    # -------------------------------------------------------------------
    # per_variant_calls.tsv
    # -------------------------------------------------------------------
    with open(args.out_variants, 'w', newline='') as fv:
        w = csv.writer(fv, delimiter='\t')
        w.writerow([
            'chrom', 'pos', 'ref', 'alt', 'gene',
            'sample_a_vaf', 'expected_ef_vaf', 'tier',
            'called_arm_a', 'called_arm_b', 'called_arm_c', 'called_arm_d',
            'called_arm_a_lofreq', 'called_arm_c_lofreq',
            'vaf_arm_a', 'vaf_arm_b', 'vaf_arm_c', 'vaf_arm_d',
            'vaf_arm_a_lofreq', 'vaf_arm_c_lofreq',
        ])
        for (chrom, pos, ref, alt), info in truth.items():
            vaf = info['vaf']
            t   = classify_tier(vaf)
            tier_totals[t] += 1
            called = {a: (chrom, pos, ref, alt) in calls[a] for a in ARMS}
            for a in ARMS:
                if called[a]:
                    tier_stats[t][a]['tp'] += 1
                else:
                    tier_stats[t][a]['fn'] += 1
            w.writerow([
                chrom, pos, ref, alt, info['gene'],
                f'{vaf:.4f}', f'{vaf / DILUTION:.4f}', t,
                *[int(called[a]) for a in ARMS],
                *[calls[a].get((chrom, pos, ref, alt), '.') for a in ARMS],
            ])

    # -------------------------------------------------------------------
    # per_arm_summary.tsv
    # -------------------------------------------------------------------
    arm_tp = {a: sum(tier_stats[t][a]['tp'] for t in [1,2,3,4]) for a in ARMS}
    arm_fn = {a: sum(tier_stats[t][a]['fn'] for t in [1,2,3,4]) for a in ARMS}

    with open(args.out_summary, 'w', newline='') as fs:
        w = csv.writer(fs, delimiter='\t')
        w.writerow([
            'arm', 'caller', 'tp', 'fn', 'fp',
            'sensitivity', 'fp_rate_per_mb', 'precision', 'f1',
            'tp_tier1', 'sens_tier1', 'tp_tier2', 'sens_tier2',
            'tp_tier3', 'sens_tier3', 'tp_tier4', 'sens_tier4',
        ])
        for a in ARMS:
            tp = arm_tp[a]; fn = arm_fn[a]; fp = fp_counts[a]
            f1, prec, sens = f1_precision_recall(tp, fp, fn)
            tier_cols = []
            for t in [1, 2, 3, 4]:
                t_tp = tier_stats[t][a]['tp']
                t_fn = tier_stats[t][a]['fn']
                t_s  = t_tp / (t_tp + t_fn) if (t_tp + t_fn) > 0 else 0.0
                tier_cols.extend([t_tp, f'{t_s:.4f}'])
            lbl, caller = arm_labels[a]
            w.writerow([
                lbl, caller, tp, fn, fp,
                f'{sens:.4f}', f'{fp / args.known_neg_mb:.2f}',
                f'{prec:.4f}', f'{f1:.4f}',
                *tier_cols,
            ])

    # -------------------------------------------------------------------
    # cross_arm_comparison.tsv — four key comparisons (Mutect2 only)
    # -------------------------------------------------------------------
    def tp_set(a):
        return frozenset(k for k in truth if k in calls[a])
    def sens(a):
        tp = arm_tp[a]; fn = arm_fn[a]
        return tp / (tp + fn) if (tp + fn) > 0 else 0.0
    def fp_mb(a):
        return fp_counts[a] / args.known_neg_mb

    comparisons = [
        ('C_vs_A', 'arm_c', 'arm_a',
         'Consensus benefit tumor-only (C vs A)'),
        ('D_vs_B', 'arm_d', 'arm_b',
         'Consensus benefit tumor-normal (D vs B)'),
        ('B_vs_A', 'arm_b', 'arm_a',
         'Paired-normal over-suppression (B vs A)'),
        ('D_vs_C', 'arm_d', 'arm_c',
         'Value of paired normal with consensus (D vs C)'),
    ]
    with open(args.out_cross_arm, 'w', newline='') as fc:
        w = csv.writer(fc, delimiter='\t')
        w.writerow([
            'comparison', 'description',
            'tp_new', 'tp_old', 'tp_delta',
            'sens_new', 'sens_old', 'sens_delta',
            'fp_mb_new', 'fp_mb_old', 'fp_mb_delta',
            'tp_gained', 'tp_lost', 'tp_shared',
        ])
        for cmp_id, a_new, a_old, desc in comparisons:
            s_new = tp_set(a_new); s_old = tp_set(a_old)
            w.writerow([
                cmp_id, desc,
                arm_tp[a_new], arm_tp[a_old], arm_tp[a_new] - arm_tp[a_old],
                f'{sens(a_new):.4f}', f'{sens(a_old):.4f}',
                f'{sens(a_new) - sens(a_old):+.4f}',
                f'{fp_mb(a_new):.2f}', f'{fp_mb(a_old):.2f}',
                f'{fp_mb(a_new) - fp_mb(a_old):+.2f}',
                len(s_new - s_old), len(s_old - s_new), len(s_new & s_old),
            ])

    # -------------------------------------------------------------------
    # fp_rate_comparison.tsv
    # -------------------------------------------------------------------
    with open(args.out_fp_rates, 'w', newline='') as ff:
        w = csv.writer(ff, delimiter='\t')
        w.writerow(['arm', 'caller', 'fp_count', 'known_neg_mb', 'fp_rate_per_mb'])
        for a in ARMS:
            lbl, caller = arm_labels[a]
            fp = fp_counts[a]
            w.writerow([lbl, caller, fp, args.known_neg_mb,
                        f'{fp / args.known_neg_mb:.2f}'])

    # Print tier summary table to stderr for visibility in logs
    print(f'\nTier summary ({len(truth)} known positives):', file=sys.stderr)
    print(f"{'Tier':<5} {'VAF (A)':<10} {'N':<5}"
          f"  {'Arm A':>12} {'Arm B':>12} {'Arm C':>12} {'Arm D':>12}",
          file=sys.stderr)
    for t in [1, 2, 3, 4]:
        n = tier_totals[t]
        vals = []
        for a in ['arm_a', 'arm_b', 'arm_c', 'arm_d']:
            tp = tier_stats[t][a]['tp']
            s  = tp / n if n > 0 else 0.0
            vals.append(f'{tp}/{n} ({s:.0%})')
        print(f"  {t:<3} {TIER_RANGES[t]:<10} {n:<5}"
              f"  {vals[0]:>12} {vals[1]:>12} {vals[2]:>12} {vals[3]:>12}",
              file=sys.stderr)

if __name__ == '__main__':
    main()
PYEOF

python3 "$PY_SCRIPT" \
  --truth-vcf       "$KNOWN_POSITIVES_VCF" \
  --arm-a           "$EVAL_DIR/arm_a_pass.tsv" \
  --arm-b           "$EVAL_DIR/arm_b_pass.tsv" \
  --arm-c           "$EVAL_DIR/arm_c_pass.tsv" \
  --arm-d           "$EVAL_DIR/arm_d_pass.tsv" \
  --arm-a-lofreq    "$EVAL_DIR/arm_a_lofreq_pass.tsv" \
  --arm-c-lofreq    "$EVAL_DIR/arm_c_lofreq_pass.tsv" \
  --fp-arm-a        "$fp_arm_a" \
  --fp-arm-b        "$fp_arm_b" \
  --fp-arm-c        "$fp_arm_c" \
  --fp-arm-d        "$fp_arm_d" \
  --fp-arm-a-lofreq "$fp_arm_a_lofreq" \
  --fp-arm-c-lofreq "$fp_arm_c_lofreq" \
  --known-neg-mb    "$KNOWN_NEG_MB" \
  --out-variants    "$EVAL_DIR/per_variant_calls.tsv" \
  --out-summary     "$EVAL_DIR/per_arm_summary.tsv" \
  --out-cross-arm   "$EVAL_DIR/cross_arm_comparison.tsv" \
  --out-fp-rates    "$EVAL_DIR/fp_rate_comparison.tsv"

rm -f "$PY_SCRIPT"

log "Script 17 complete. Results in $EVAL_DIR/"
log "  per_variant_calls.tsv   — $(tail -n +2 "$EVAL_DIR/per_variant_calls.tsv" | wc -l) variants"
log "  per_arm_summary.tsv"
log "  cross_arm_comparison.tsv"
log "  fp_rate_comparison.tsv"
