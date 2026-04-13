# ctDNA Variant Calling Pipeline Benchmark

A reproducible benchmarking pipeline for somatic variant detection in 
circulating tumour DNA (ctDNA), comparing coordinate-based consensus 
calling against standard duplicate-marking approaches across four 
experimental arms.

Built on publicly available SEQC2 reference data with a published 
ground truth, this pipeline demonstrates end-to-end ctDNA variant 
calling from raw FASTQ through filtered VCF, with quantified 
sensitivity and false positive rate across clinically relevant 
VAF tiers.

---

## Background and Motivation

### Why ctDNA variant calling is hard

Circulating tumour DNA represents a small and variable fraction of 
total cell-free DNA in plasma — typically 0.1–5% in early-stage or 
MRD monitoring contexts. Detecting somatic mutations at these allele 
frequencies requires sequencing at extreme depth (>10,000x) and 
rigorous error suppression, because the sequencing error rate 
(~0.1–0.5% per base) is comparable to or greater than the variant 
allele frequency being measured.

### The role of UMIs

Unique Molecular Identifiers (UMIs) are short random oligonucleotide 
sequences ligated to each DNA fragment before PCR amplification. 
Because every original molecule receives a unique tag, reads sharing 
the same UMI originate from the same molecule. This allows:

- **PCR duplicate identification** without relying on mapping 
  coordinates alone
- **Consensus base calling** across multiple reads from the same 
  molecule, suppressing both PCR errors and sequencing errors
- **Accurate duplicate rate estimation** even when two distinct 
  molecules happen to map to identical coordinates

At sub-1% VAF, the error rate reduction from UMI consensus is the 
difference between a detectable signal and noise.

### The UMI archival problem

UMI sequences are frequently stripped during data archival. NCBI SRA's 
ETL normalisation pipeline discards UMI adapter sequences when labs 
submit pre-processed FASTQs rather than raw BCL or unprocessed BAM 
files. This was confirmed for the SEQC2 BRP dataset used here via 
`vdb-dump` (FMT: FASTQ, LDR: latf-load) and ENA filereport 
(empty submitted_ftp field).

This pipeline addresses UMI loss by reconstructing molecule identity 
from fragment coordinates — a valid proxy for cfDNA data where 
nucleosomal protection creates defined fragment endpoints — and 
quantifies the downstream impact on variant calling performance.

---

## Experimental Design

### Data

| Sample | SRR | Description | Tumor Fraction |
|--------|-----|-------------|----------------|
| Df (tumor) | SRR13200966 | SEQC2 BRP v4, Site 25, 25ng, LIB1 | ~20% (1:4 dilution) |
| Bf (normal) | SRR13201012 | SEQC2 BRP v4, Site 25, 25ng, LIB1 | 0% |

Both samples were sequenced with the Burning Rock Lung Plasma v4 
ctDNA assay targeting ~70 clinically relevant oncology genes across 
~53kb of coding sequence. Df contains Sample A cancer cell line 
variants diluted 1:4 into Sample B normal background.

### Four pipeline arms

| Arm | Duplicate handling | Calling mode | Callers |
|-----|-------------------|--------------|---------|
| A | Picard MarkDuplicates | Tumor-only | Mutect2, LoFreq |
| B | Picard MarkDuplicates | Tumor-normal | Mutect2 |
| C | fgbio coordinate consensus | Tumor-only | Mutect2, LoFreq |
| D | fgbio coordinate consensus | Tumor-normal | Mutect2 |

### Truth set

The SEQC2 known positives VCF (42,317 variants) was intersected with 
the BRP panel (hg38 liftover coordinates) and filtered to remove 
germline variants present in gnomAD at AF ≥ 0.001. This yielded 
**86 somatic candidates** — the evaluable truth set.

| Category | Count | Notes |
|----------|-------|-------|
| Common germline (gnomAD >1%) | 124 | Correctly filtered by Mutect2 |
| Rare germline (gnomAD 0.001–1%) | 18 | Partially filtered |
| Somatic candidates | 86 | Evaluation truth set |

All 86 somatic candidates fall in Tier 3 (Sample A VAF 5–20%, 
expected ~1–4% in Df) or Tier 4 (Sample A VAF <5%, expected <1% 
in Df).

---

## Key Results

### Sensitivity by pipeline arm

| Arm | Overall | Tier 3 (1–4% VAF) | Tier 4 (<1% VAF) | FP/Mb |
|-----|---------|-------------------|------------------|-------|
| A — Standard tumor-only | 17.4% | 39.3% | 6.9% | 132 |
| B — Standard tumor-normal | 2.3% | 7.1% | 0.0% | 0 |
| C — Consensus tumor-only | 69.8% | 71.4% | 69.0% | 470 |
| D — Consensus tumor-normal | 52.3% | 32.1% | 62.1% | 19 |

### Headline findings

**1. Consensus calling dramatically improves low-VAF sensitivity.**
Coordinate-based consensus (Arm C) recovers 70% of somatic candidates 
vs 17% for standard calling (Arm A) — a 4x improvement driven 
entirely by the 3.3x reduction in per-base error rate achieved through 
family-level consensus (0.00049 → 0.00016).

**2. Paired-normal over-suppression is severe without consensus.**
Standard tumor-normal calling (Arm B) achieves near-perfect specificity 
(0 FP/Mb) but at the cost of 85% sensitivity loss relative to 
tumor-only (Arm A). Both samples share Sample B germline background, 
causing Mutect2 to suppress real somatic variants present in the 
matched normal. This is the clinical CHIP problem quantified empirically.

**3. Consensus + paired normal achieves the best clinical balance.**
Arm D recovers 52% of somatic candidates with only 1 false positive 
(19 FP/Mb, F1=0.68, precision=0.978). Compared to standard 
tumor-normal (Arm B, F1=0.046), consensus error reduction rescues 
43 additional true positives while adding only 1 FP — demonstrating 
that coordinate-based consensus is sufficient to make low-VAF 
tumor-normal calling clinically viable.

**4. LoFreq is highly sensitive but requires consensus for specificity.**
LoFreq achieves 93% sensitivity in both Arms A and C but at 4005 and 
2388 FP/Mb respectively. Consensus reduces LoFreq's FP burden by 40% 
while maintaining sensitivity — suggesting a combined Mutect2+LoFreq 
approach with consensus as a pre-filter may be optimal for MRD 
detection.

**5. UMI preservation is the critical unresolved limitation.**
True UMI-aware calling would further reduce false positives by 
resolving the rare case of two distinct molecules with identical 
fragment coordinates — a limitation of the coordinate-consensus 
approach used here. The performance gap between Arms C/D and a 
true UMI-aware pipeline represents the quantifiable cost of UMI 
stripping during SRA archival.

---

## Repository Structure

```
├── scripts/           # Pipeline scripts 00–18 + run_all.sh
├── ref/               # Reference files and truth sets
│   ├── BRP2/          # BRP panel BED (hg19 and hg38 liftover)
│   └── SampleA_ref/   # SEQC2 truth VCFs and known negatives
├── docs/              # Design notes and methodology
├── results/           # Pipeline outputs (VCFs, metrics, REPORT.md)
├── environment.yml    # Conda environment specification
└── CLAUDE.md          # Claude Code session context
```

---

## Quickstart

```bash
# Clone and create environment
git clone https://github.com/gavinoliver-prog/umi-aware-seqc2-ctdna-variant-calling
cd umi-aware-seqc2-ctdna-variant-calling
conda env create -f environment.yml
conda activate ctdna2

# Run setup (downloads fgbio 2.3.0 JAR, verifies tools)
bash scripts/01_setup_env.sh

# Download data (FASTQs ~22GB total)
bash scripts/02_download_data.sh

# Run full pipeline (~24 hours on 8-core instance)
bash scripts/run_all.sh
```

Full results are in `results/REPORT.md`.

---

## Technical Notes

**Why coordinate-based consensus rather than true UMI calling:**
UMIs were stripped during SRA submission of the SEQC2 BRP dataset. 
Fragment coordinate grouping is used as a valid proxy for cfDNA 
(nucleosomal endpoints are biologically defined) but cannot 
distinguish two distinct molecules with identical endpoints. See 
`docs/design_notes.md` for the full investigation.

**Why Df rather than Ef:**
Initial analysis with Ef (1:24 dilution, ~4% tumor fraction) showed 
all evaluable somatic variants at <0.8% expected VAF — below the 
reliable detection threshold for standard calling. Df (1:4 dilution, 
~20% tumor fraction) places Tier 3 variants at 1–4% VAF, within the 
detectable range. See `docs/design_notes.md`.

**Why the truth set was reduced from 228 to 86 variants:**
54% of the original SEQC2 known positives overlapping the BRP panel 
are common germline SNPs (gnomAD AF >1%) correctly filtered by 
Mutect2's germline model. Evaluating against these would penalise 
correct germline filtering. The refined somatic truth set excludes 
variants present in gnomAD at AF ≥ 0.001.

---

## Full Report

See [results/REPORT.md](results/REPORT.md) for the complete 
benchmark report including consensus metrics, per-tier sensitivity 
tables, cross-arm comparisons, and limitations.

---

*Pipeline developed as a portfolio project demonstrating ctDNA 
bioinformatics methodology. SEQC2 data is publicly available at 
PRJNA715852.*
