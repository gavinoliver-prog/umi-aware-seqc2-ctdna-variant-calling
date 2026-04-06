#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Generate results/REPORT.md — a human-readable summary of all pipeline results.
#
# Reads TSV outputs from scripts 15 (consensus metrics) and 17 (evaluation).
# Python 3 stdlib only — no bioinformatics tools required.

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_cmd python3

METRICS_DIR="$RESULTS_DIR/consensus_metrics"
EVAL_DIR="$RESULTS_DIR/evaluation"
REPORT="$RESULTS_DIR/REPORT.md"

require_file "$METRICS_DIR/summary_metrics.tsv"
require_file "$METRICS_DIR/family_size_distribution.tsv"
require_file "$METRICS_DIR/error_rate_comparison.tsv"
require_file "$EVAL_DIR/per_arm_summary.tsv"
require_file "$EVAL_DIR/per_variant_calls.tsv"
require_file "$EVAL_DIR/cross_arm_comparison.tsv"
require_file "$EVAL_DIR/fp_rate_comparison.tsv"

if [[ -f "$REPORT" ]]; then
  log "Report exists: $REPORT — skipping."
  exit 0
fi

log "Generating report..."

python3 - \
  --metrics-dir "$METRICS_DIR" \
  --eval-dir    "$EVAL_DIR" \
  --output      "$REPORT" \
  --min-reads   "$UMI_MIN_READS" \
  --min-bq      "$UMI_MIN_BASE_QUALITY" \
  --tumor-run   "$TUMOR_RUN" \
  --normal-run  "$NORMAL_RUN" \
  <<'PYEOF'
import argparse, csv, sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--metrics-dir', required=True)
ap.add_argument('--eval-dir',    required=True)
ap.add_argument('--output',      required=True)
ap.add_argument('--min-reads',   type=int, default=2)
ap.add_argument('--min-bq',      type=int, default=20)
ap.add_argument('--tumor-run',   default='SRR13200999')
ap.add_argument('--normal-run',  default='SRR13201012')
args = ap.parse_args()

mdir = Path(args.metrics_dir)
edir = Path(args.eval_dir)
out  = Path(args.output)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def read_tsv(path):
    """Return list of dicts from a TSV file."""
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def md_table(rows, headers=None):
    """Render a list-of-dicts (or list-of-lists) as a Markdown table."""
    if not rows:
        return '_No data._\n'
    if headers is None:
        headers = list(rows[0].keys()) if isinstance(rows[0], dict) else None
    if isinstance(rows[0], dict):
        data = [[str(r.get(h, '')) for h in headers] for r in rows]
    else:
        data = [[str(c) for c in r] for r in rows]
    col_widths = [max(len(h), max((len(r[i]) for r in data), default=0))
                  for i, h in enumerate(headers)]
    def fmt_row(cells):
        return '| ' + ' | '.join(c.ljust(w) for c, w in zip(cells, col_widths)) + ' |'
    sep = '| ' + ' | '.join('-' * w for w in col_widths) + ' |'
    lines = [fmt_row(headers), sep] + [fmt_row(r) for r in data]
    return '\n'.join(lines) + '\n'

def fmt_rate(val, decimals=6):
    try:
        return f'{float(val):.{decimals}f}'
    except (ValueError, TypeError):
        return str(val)

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

summary_rows = read_tsv(mdir / 'summary_metrics.tsv')
family_df    = read_tsv(mdir / 'family_size_distribution.tsv')
error_df     = read_tsv(mdir / 'error_rate_comparison.tsv')

per_arm   = read_tsv(edir / 'per_arm_summary.tsv')
per_var   = read_tsv(edir / 'per_variant_calls.tsv')
cross_arm = read_tsv(edir / 'cross_arm_comparison.tsv')
fp_rate   = read_tsv(edir / 'fp_rate_comparison.tsv')

# -------------------------------------------------------------------------
# Lookup helpers matched to actual TSV schemas
# -------------------------------------------------------------------------

def get_metric(metric, col):
    """summary_metrics.tsv: columns are metric, tumor, normal."""
    for r in summary_rows:
        if r.get('metric') == metric:
            return r.get(col, 'N/A')
    return 'N/A'

def arm_sens(arm_label, caller='Mutect2'):
    """per_arm_summary.tsv: arm in {A,B,C,D}, caller in {Mutect2,LoFreq}."""
    for r in per_arm:
        if r.get('arm') == arm_label and r.get('caller') == caller:
            return r.get('sensitivity', 'N/A')
    return 'N/A'

def get_cross(comparison, field):
    """cross_arm_comparison.tsv: fields include sens_new, sens_old, tp_gained, etc."""
    for r in cross_arm:
        if r.get('comparison') == comparison:
            return r.get(field, 'N/A')
    return 'N/A'

def arm_fp_rate(arm_label, caller='Mutect2'):
    """fp_rate_comparison.tsv: arm in {A,B,C,D}, column is fp_rate_per_mb."""
    for r in fp_rate:
        if r.get('arm') == arm_label and r.get('caller') == caller:
            return r.get('fp_rate_per_mb', 'N/A')
    return 'N/A'

# Derived summary numbers
t_pre_reads  = get_metric('pre_consensus_total_reads',  'tumor')
t_dup        = get_metric('markdup_duplicate_rate',     'tumor')
t_post_reads = get_metric('post_consensus_total_reads', 'tumor')
t_comp       = get_metric('compression_ratio',          'tumor')
t_pre_err    = get_metric('error_rate_pre_consensus',   'tumor')
t_post_err   = get_metric('error_rate_post_consensus',  'tumor')

n_pre_reads  = get_metric('pre_consensus_total_reads',  'normal')
n_dup        = get_metric('markdup_duplicate_rate',     'normal')
n_post_reads = get_metric('post_consensus_total_reads', 'normal')
n_comp       = get_metric('compression_ratio',          'normal')
n_pre_err    = get_metric('error_rate_pre_consensus',   'normal')
n_post_err   = get_metric('error_rate_post_consensus',  'normal')

# -------------------------------------------------------------------------
# Build report
# -------------------------------------------------------------------------

lines = []
L = lines.append

def section(title, level=2):
    L('')
    L('#' * level + ' ' + title)
    L('')

# =========================================================================
# 1. Title and executive summary
# =========================================================================
L('# ctDNA Variant Calling Pipeline Benchmark')
L('')
L('**Coordinate-based consensus calling vs standard duplicate marking**  ')
L('**for low-VAF somatic variant detection in ctDNA**')
L('')
L('---')

section('Executive Summary')
L('This benchmark evaluates four variant calling strategies on SEQC2 Burning Rock ')
L('Lung Plasma (BRP v4) reference material. The tumor sample (Ef) is a 1:24 dilution ')
L('of cancer cell line into normal background, yielding an expected tumor fraction of ')
L('~4%. Against 228 panel-overlapping known positive variants and 53,185 bp of known ')
L('negative space, the four arms reveal:')
L('')

a_sens = arm_sens('A'); b_sens = arm_sens('B')
c_sens = arm_sens('C'); d_sens = arm_sens('D')
a_fp   = arm_fp_rate('A'); c_fp = arm_fp_rate('C')

L(f'- **Arm A (standard tumor-only):** {a_sens} overall sensitivity, {a_fp} FP/Mb')
L(f'- **Arm B (standard tumor-normal):** {b_sens} overall sensitivity — paired-normal '
  f'over-suppression removes variants shared with normal background')
L(f'- **Arm C (consensus tumor-only):** {c_sens} overall sensitivity, {c_fp} FP/Mb')
L(f'- **Arm D (consensus tumor-normal):** {d_sens} overall sensitivity')
L('')
L('Key findings are described in sections 6–8.')

# =========================================================================
# 2. Dataset and experimental design
# =========================================================================
section('Dataset and Experimental Design')
L(f'**Tumor:** `{args.tumor_run}` (Ef — 1:24 Sample A cancer cell line : Sample B normal)  ')
L(f'**Normal:** `{args.normal_run}` (Bf — pure Sample B normal)  ')
L('**Panel:** Burning Rock BRP v4 (~53 kb coding regions, hg38)')
L('')
L('### Four Pipeline Arms')
L('')
arm_table = [
    {'Arm': 'A', 'Strategy': 'Standard tumor-only',   'Callers': 'Mutect2 + LoFreq'},
    {'Arm': 'B', 'Strategy': 'Standard tumor-normal',  'Callers': 'Mutect2'},
    {'Arm': 'C', 'Strategy': 'Consensus tumor-only',   'Callers': 'Mutect2 + LoFreq'},
    {'Arm': 'D', 'Strategy': 'Consensus tumor-normal', 'Callers': 'Mutect2'},
]
L(md_table(arm_table))
L('Arms A and B share the same tumor BAM (duplicate-marked). '
  'Arms C and D share the same consensus tumor BAM.')

# =========================================================================
# 3. Methodological note: coordinate consensus vs true UMI calling
# =========================================================================
section('Methodological Note: Coordinate Consensus vs True UMI Calling')
L('UMIs were stripped during SRA submission of the SEQC2 BRP dataset '
  '(confirmed via `vdb-dump`). Arms C and D therefore implement '
  '**coordinate-based consensus calling**, not true UMI-aware calling.')
L('')
L('| Approach | Duplicate identification | Errors corrected |')
L('|---|---|---|')
L('| True UMI calling | UMI sequence ligated before PCR | PCR duplicates + sequencing error |')
L('| Coordinate consensus (this pipeline) | chr + start + end + strand | Sequencing error only |')
L('')
L('`fgbio GroupReadsByUmi --strategy identity` is used because synthetic UMIs are derived '
  'from fragment coordinates — adjacency/edit-distance correction would be meaningless '
  '(a 1-base UMI difference encodes a different coordinate group, not a PCR error).')
L('')
L('This limitation is itself a finding: it quantifies the cost of UMI data loss during '
  'archival and motivates proper UMI preservation in future studies.')

# =========================================================================
# 4. Consensus calling metrics
# =========================================================================
section('Consensus Calling Metrics')

L('### Read counts and compression')
L('')
metrics_table = [
    {'Metric': 'Pre-consensus reads',       'Tumor': t_pre_reads,  'Normal': n_pre_reads},
    {'Metric': 'Duplicate rate',            'Tumor': t_dup,         'Normal': n_dup},
    {'Metric': 'Post-consensus reads',      'Tumor': t_post_reads,  'Normal': n_post_reads},
    {'Metric': 'Compression ratio',         'Tumor': t_comp,        'Normal': n_comp},
    {'Metric': 'Error rate pre-consensus',  'Tumor': fmt_rate(t_pre_err),  'Normal': fmt_rate(n_pre_err)},
    {'Metric': 'Error rate post-consensus', 'Tumor': fmt_rate(t_post_err), 'Normal': fmt_rate(n_post_err)},
]
L(md_table(metrics_table))
L('')
L('_Error rate proxy: non-reference allele rate at known-negative loci (BQ ≥ 20). '
  'Absolute values are inflated by Sample B germline variants present in both samples; '
  'the pre- vs post-consensus comparison is valid._')
L('')

# Family size distribution — show up to 15 rows per sample
L('### Family size distribution (fgbio GroupReadsByUmi)')
L('')
L('Family size = number of read pairs sharing identical fragment coordinates '
  '(the "duplicate group" size). Families of size 1 are singletons; '
  f'`--min-reads {args.min_reads}` discards these in consensus calling.')
L('')
for sample in ('tumor', 'normal'):
    sample_rows = [r for r in family_df if r.get('sample') == sample][:15]
    if sample_rows:
        L(f'**{sample.capitalize()}**')
        L('')
        L(md_table(sample_rows, headers=['family_size', 'count', 'fraction']))
        L('')

# =========================================================================
# 5. Variant calling results
# =========================================================================
section('Variant Calling Results')

L('### Per-arm summary')
L('')
# per_arm_summary has one row per arm+caller — show key columns
if per_arm:
    display_cols = ['arm', 'caller', 'tp', 'fn', 'fp',
                    'sensitivity', 'fp_rate_per_mb', 'precision', 'f1']
    L(md_table(per_arm, headers=display_cols))
L('')

L('### Sensitivity by VAF tier')
L('')
L('Expected VAF in Ef = Sample A VAF / 25 (1:24 dilution)')
L('')
# Tier data is stored as flat columns in per_arm_summary.tsv:
# tp_tier1, sens_tier1, tp_tier2, sens_tier2, tp_tier3, sens_tier3, tp_tier4, sens_tier4
tier_label = {1: '>=50% (>2% in Ef)', 2: '20-50% (0.8-2%)',
              3: '5-20% (0.2-0.8%)',  4: '<5% (<0.2%)'}
tier_rows = []
for r in per_arm:
    arm = r.get('arm', '?'); caller = r.get('caller', '?')
    for t in [1, 2, 3, 4]:
        tier_rows.append({
            'arm':         arm,
            'caller':      caller,
            'tier':        t,
            'vaf_range':   tier_label[t],
            'tp':          r.get(f'tp_tier{t}', 'N/A'),
            'sensitivity': r.get(f'sens_tier{t}', 'N/A'),
        })
if tier_rows:
    L(md_table(tier_rows))
L('')

# =========================================================================
# 6. Key finding 1: consensus benefit
# =========================================================================
section('Key Finding 1: Consensus Calling Benefit (Arm C vs A)')
c_vs_a_gained = get_cross('C_vs_A', 'tp_gained')
c_vs_a_lost   = get_cross('C_vs_A', 'tp_lost')
c_vs_a_shared = get_cross('C_vs_A', 'tp_shared')
c_vs_a_sens_a = get_cross('C_vs_A', 'sens_old')
c_vs_a_sens_c = get_cross('C_vs_A', 'sens_new')
L('Consensus calling (Arm C) vs standard duplicate marking (Arm A) in the tumor-only setting:')
L('')
L(f'- **TP shared:** {c_vs_a_shared} variants called by both arms')
L(f'- **TP gained by C:** {c_vs_a_gained} variants called only by Arm C')
L(f'- **TP lost by C:** {c_vs_a_lost} variants called only by Arm A')
L(f'- **Overall sensitivity:** Arm A = {c_vs_a_sens_a}, Arm C = {c_vs_a_sens_c}')
L('')
L(f'FP rate: Arm A = {arm_fp_rate("A")} FP/Mb, Arm C = {arm_fp_rate("C")} FP/Mb (Mutect2)')
L('')
L('Consensus base quality improvement enables Mutect2 and LoFreq to better distinguish '
  'true low-VAF variants from sequencing noise at the cost of reduced read depth '
  'after family collapsing.')

# =========================================================================
# 7. Key finding 2: paired-normal over-suppression
# =========================================================================
section('Key Finding 2: Paired-Normal Over-Suppression (Arm B vs A)')
b_vs_a_gained = get_cross('B_vs_A', 'tp_gained')
b_vs_a_lost   = get_cross('B_vs_A', 'tp_lost')
b_vs_a_sens_a = get_cross('B_vs_A', 'sens_old')
b_vs_a_sens_b = get_cross('B_vs_A', 'sens_new')
L('Both tumor (Ef) and normal (Bf) samples contain Sample B germline background. '
  'Mutect2 tumor-normal calling suppresses variants present in the normal — which '
  'includes real tumor variants that happen to overlap the shared germline signal.')
L('')
L(f'- **Overall sensitivity:** Arm A = {b_vs_a_sens_a}, Arm B = {b_vs_a_sens_b}')
L(f'- **TP lost vs Arm A:** {b_vs_a_lost} variants suppressed by the paired normal')
L(f'- **TP gained vs Arm A:** {b_vs_a_gained} variants rescued from FP filtering by the normal')
L('')
L('This over-suppression is a known clinical problem related to CHIP and shared '
  'germline signal. Quantifying it empirically is a key result of this pipeline.')

# =========================================================================
# 8. Key finding 3: FP rate comparison
# =========================================================================
section('Key Finding 3: False Positive Rate Comparison')
L(md_table(fp_rate))
L('')
L('FP rate is computed over 53,185 bp of known-negative space. '
  'LoFreq FP rates include all calls (LoFreq does not use a PASS filter by default); '
  'Mutect2 rates include PASS-only calls after FilterMutectCalls.')

# =========================================================================
# 9. Limitations and future work
# =========================================================================
section('Limitations and Future Work')
L('1. **No true UMIs:** Coordinate-based grouping assumes identical-endpoint '
  'read pairs are PCR duplicates. Two distinct molecules with coincidentally '
  'identical endpoints are conflated — this is rare but real, and cannot be '
  'disambiguated without the original UMI sequence.')
L('')
L('2. **Low tumor fraction:** At ~4% tumor fraction and 1:24 dilution, most Tier 3 '
  'and Tier 4 variants are below 1% VAF. Higher-coverage sequencing or higher '
  'tumor fraction samples would be needed to assess performance in the >2% VAF range.')
L('')
L('3. **Germline contamination in error rate proxy:** Known-negative error rates '
  'are inflated by Sample B germline variants. True sequencing error rates are lower; '
  'a germline-filtered analysis would provide cleaner baseline estimates.')
L('')
L('4. **Single replicate:** Performance estimates would be more robust with '
  'multiple Ef replicates at different tumor fractions (BRP provides several).')
L('')
L('5. **UMI archival:** The primary recommendation from this study is that UMI '
  'sequences should be preserved in submitted FASTQ files. The SRA/ENA submission '
  'pipeline should be configured to retain RX tags or UMI-prefix read names.')

# =========================================================================
# Write
# =========================================================================
report_text = '\n'.join(lines) + '\n'
out.write_text(report_text)
print(f'Report written: {out}', file=sys.stderr)
PYEOF

log "Report written: $REPORT"
log "Script 18 complete."
