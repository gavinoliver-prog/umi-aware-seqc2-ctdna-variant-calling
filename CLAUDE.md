> This file provides context for Claude Code (Anthropic's agentic 
> coding assistant) during development sessions. It documents the 
> scientific design decisions, methodological rationale, and coding 
> standards for this pipeline. It is preserved in the repository as 
> supplementary documentation of the development process.

# ctDNA Variant Calling Pipeline — Claude Code Project Context

## Project Summary

Benchmarking coordinate-based consensus calling vs standard duplicate marking
for low-VAF somatic variant detection in ctDNA data. Uses SEQC2 BRP 
(Burning Rock Lung Plasma v4) reference samples with published ground truth.

This is a bioinformatics portfolio project demonstrating:
1. The limitations of standard duplicate marking at low VAF
2. The benefit of consensus calling (fgbio) for error reduction
3. The paired-normal over-suppression effect in ctDNA tumor-normal calling
4. Why UMI preservation during data archival matters

## Critical Methodological Note

**UMIs were stripped during SRA submission of this dataset** (confirmed via
vdb-dump showing FMT: FASTQ and ENA filereport showing empty submitted_ftp).

Therefore Arms C and D implement **coordinate-based consensus calling**, NOT
true UMI-aware calling. The distinction is important:

- True UMI calling: groups reads by UMI sequence ligated before PCR
  Proves common origin, resolves identical-coordinate ambiguity
  
- Coordinate consensus (this pipeline): groups reads by chr+start+end+strand
  Same assumption as Picard MarkDuplicates for duplicate identification
  Adds value through consensus base quality improvement over read discarding
  Cannot resolve two distinct molecules with identical endpoints (rare but real)

fgbio GroupReadsByUmi is used with --strategy identity because synthetic UMIs
are derived from coordinates — adjacency/edit-distance correction would be
redundant and misleading.

This limitation is itself a finding: it quantifies the cost of UMI data loss
during archival and motivates proper UMI preservation in future studies.

umi-tools is NOT used — it would be redundant with fgbio for coordinate-based
grouping since UMIs are derived from the same coordinates used for grouping.

## Experimental Design

### Samples
- Tumor:  SRR13200999 (Ef — 1:24 tumor:normal, ~4% tumor fraction)
- Normal: SRR13201012 (Bf — 0:1 tumor:normal, pure normal)

Both are SEQC2 BRP replicates. Ef contains Sample A cancer cell line variants
diluted 1:24 into Sample B normal background.

Expected VAF in Ef = Sample A VAF / 25

### Four Pipeline Arms

Arm A: Standard tumor-only
       Ef → BWA-MEM → Picard MarkDuplicates → Mutect2
       Baseline: standard of care, no error correction

Arm B: Standard tumor-normal
       Ef+Bf → BWA-MEM → Picard MarkDuplicates → Mutect2
       Expected: paired-normal over-suppression because Ef and Bf share
       Sample B germline background — real tumor variants will be suppressed.
       This over-suppression is a known clinical problem (related to CHIP
       and shared germline signal) and observing it empirically is a key result.

Arm C: Consensus tumor-only
       Ef → BWA-MEM → fgbio coordinate consensus → Mutect2 + LoFreq
       Tests whether consensus base quality improvement aids sensitivity
       without the complication of paired-normal suppression

Arm D: Consensus tumor-normal
       Ef+Bf → BWA-MEM → fgbio coordinate consensus → Mutect2
       Tests whether consensus calling partially rescues paired-normal
       over-suppression by reducing noise floor

### Key Comparisons
- Arm C vs A: benefit of consensus calling in tumor-only setting
- Arm D vs B: benefit of consensus calling in tumor-normal setting
- Arm B vs A: paired-normal over-suppression (variants lost, quantified)
- Arm D vs C: added value of paired normal when using consensus calling

## Truth Set

### Known Positives
228 variants overlapping BRP panel (hg38 liftover coordinates)
File: ref/SampleA_ref/KnownPositives_hg38.on_target.vcf.gz
VAF field in INFO = Sample A VAF (before dilution)

VAF stratification:
  Tier 1: Sample A VAF >=50%  → >2% expected in Ef    (30 variants)
  Tier 2: Sample A VAF 20-50% → 0.8-2% in Ef          (39 variants)
  Tier 3: Sample A VAF 5-20%  → 0.2-0.8% in Ef        (76 variants)
  Tier 4: Sample A VAF <5%    → <0.2% in Ef            (83 variants)

Key genes: EGFR, TP53, BRCA1, BRCA2, KRAS, BRAF, ALK, PIK3CG, KEAP1,
           STK11, CDKN2A, APC, RET, ERBB2, PTEN, MET, FGFR3, POLE

### Known Negatives
53,185 bp overlapping BRP panel
File: ref/SampleA_ref/KnownNegatives_hg38.bed.gz
FP rate = FP calls / 0.053185 Mb (variants per Mb of known negative space)

## Directory Structure

PROJECT_ROOT=/home/ubuntu/projects/ctdna2
DATA_DIR=$PROJECT_ROOT/data
FASTQ_DIR=$DATA_DIR/fastq          # FASTQs already downloaded
REF_DIR=$PROJECT_ROOT/ref
RESULTS_DIR=$PROJECT_ROOT/results
BAM_DIR=$RESULTS_DIR/bam
VCF_DIR=$RESULTS_DIR/vcf
TMP_DIR=$DATA_DIR/tmp
ARM_A_DIR=$RESULTS_DIR/arm_a_standard_tumor_only
ARM_B_DIR=$RESULTS_DIR/arm_b_standard_tumor_normal
ARM_C_DIR=$RESULTS_DIR/arm_c_consensus_tumor_only
ARM_D_DIR=$RESULTS_DIR/arm_d_consensus_tumor_normal

## Key Reference Files

REF_FASTA:        ref/Homo_sapiens_assembly38.fasta (BWA index complete)
TARGET_BED:       ref/BRP2/BRP2_liftover_hg19tohg38.bed
KNOWN_POS_VCF:    ref/SampleA_ref/KnownPositives_hg38.on_target.vcf.gz
KNOWN_NEG_BED:    ref/SampleA_ref/KnownNegatives_hg38.bed.gz
gnomAD:           ref/af-only-gnomad.hg38.vcf.gz (may not exist — skip gracefully)
PON:              ref/1000g_pon.hg38.vcf.gz (may not exist — skip gracefully)

## Conda Environment

Name: ctdna2
Installed: bwa, samtools, gatk4, picard, bcftools, fastqc, multiqc
MISSING — must install: fgbio, lofreq

## Existing Scripts (scripts/ directory)

00_config.sh          Central config — all scripts source this first
01_setup_env.sh       Environment setup
02_download_data.sh   Data download — FASTQs already present
03_subsample_fastq.sh Subsampling
04_prepare_reference.sh Reference prep — BWA index already complete
05_fastqc_multiqc.sh  QC
06_align_and_qc.sh    Alignment — written for TST170, needs BRP adaptation
07_call_mutect_tumor_only.sh    Tumor-only calling
08_coverage_and_truth.sh        Coverage analysis
09_call_mutect_tumor_normal.sh  Tumor-normal calling
10_compare_results.sh           Comparison — needs rewrite for 4-arm design
run_all.sh            Orchestration

## Required Changes to Existing Scripts

### 00_config.sh — add these variables:
TUMOR_RUN=SRR13200999         # Ef (was MIX124_RUN)
NORMAL_RUN=SRR13201012        # Bf (was MIX01_RUN)
BRP_DIR=$REF_DIR/BRP2
TARGET_BED=$BRP_DIR/BRP2_liftover_hg19tohg38.bed
KNOWN_POSITIVES_VCF=$REF_DIR/SampleA_ref/KnownPositives_hg38.on_target.vcf.gz
KNOWN_NEGATIVES_BED=$REF_DIR/SampleA_ref/KnownNegatives_hg38.bed.gz
ON_TARGET_TRUTH_VCF=$REF_DIR/SampleA_ref/KnownPositives_hg38.on_target.vcf.gz
ARM_A_DIR=$RESULTS_DIR/arm_a_standard_tumor_only
ARM_B_DIR=$RESULTS_DIR/arm_b_standard_tumor_normal
ARM_C_DIR=$RESULTS_DIR/arm_c_consensus_tumor_only
ARM_D_DIR=$RESULTS_DIR/arm_d_consensus_tumor_normal
UMI_MIN_READS=2
UMI_MIN_BASE_QUALITY=20
FGBIO_JAR=""   # set after install in 11_install_missing_tools.sh

### 06_align_and_qc.sh:
Update TARGET_BED reference to use BRP panel not TST170

## New Scripts to Create

### 11_install_missing_tools.sh
Install fgbio and lofreq into ctdna2 conda environment.
Idempotent — check if already installed first.
Set FGBIO_JAR path and verify it works.

### 12_group_fragments_by_coordinate.sh
METHODOLOGICALLY CRITICAL — comment extensively.

This script implements coordinate-based fragment grouping as a proxy for
UMI-based duplicate identification. This is NOT equivalent to true UMI
calling — see methodological note above.

Steps:
1. Align Ef FASTQs to hg38 with BWA-MEM (temporary alignment for coord extraction)
2. Extract on-target read pairs overlapping BRP panel BED
3. Group read pairs by fragment coordinates (chr + leftmost + rightmost + strand)
   Same coordinates = assumed same original molecule (same as MarkDuplicates)
4. Assign deterministic synthetic 8bp UMI per coordinate group
   Use MD5 hash of "chr:start:end:strand" truncated to 8 ACGT bases
   Deterministic = reproducible across runs
   These UMIs encode NO additional information beyond the coordinates
5. Write UMI-prefixed FASTQs (read structure: 8M+T +T) for fgbio pipeline
6. Write standard FASTQs (no UMI prefix) for Arms A/B — same reads, fair comparison
7. Output QC: read pairs processed, on-target rate, unique fragments,
   family size distribution, estimated duplicate rate
8. Repeat for Bf (normal) sample

### 13_align_arms_ab.sh
Standard pipeline — Arms A and B.
Input: standard FASTQs from script 12
BWA-MEM → Picard MarkDuplicates → collect metrics
Arms A and B share same tumor BAM — calling differs not alignment
Output: ARM_A_DIR/tumor.markdup.bam, ARM_B_DIR/normal.markdup.bam

### 14_align_arms_cd.sh
Consensus pipeline — Arms C and D.
Input: UMI-prefixed FASTQs from script 12

fgbio pipeline:
  FastqToBam (read structure 8M+T +T — extract synthetic UMI to RX tag)
  → BWA-MEM via: samtools fastq | bwa mem | fgbio ZipperBams
  → fgbio GroupReadsByUmi --strategy identity
    (identity not adjacency — UMIs are coordinate-derived, correction meaningless)
  → fgbio CallMolecularConsensusReads --min-reads $UMI_MIN_READS
  → fgbio FilterConsensusReads --min-base-quality $UMI_MIN_BASE_QUALITY
  → re-align consensus reads with BWA-MEM

Arms C and D share same consensus tumor BAM
Output: ARM_C_DIR/tumor.consensus.bam, ARM_D_DIR/normal.consensus.bam

### 15_collect_consensus_metrics.sh
Quantify the effect of consensus calling:
- Pre-consensus: total reads, duplicate rate
- Post-consensus: read count, compression ratio (pre/post)
- Family size distribution from fgbio GroupReadsByUmi metrics
- Error rate proxy: variant rate in known-negative regions before/after consensus
  (bcftools mpileup on known negatives BED, count non-reference bases)
Output: results/consensus_metrics/family_size_distribution.tsv
        results/consensus_metrics/summary_metrics.tsv
        results/consensus_metrics/error_rate_comparison.tsv

### 16_call_variants_all_arms.sh
All four arms, restricted to TARGET_BED with -L flag throughout.

Arm A: Mutect2 tumor-only on ARM_A_DIR/tumor.markdup.bam
       LoFreq tumor-only on same BAM
Arm B: Mutect2 tumor-normal
       tumor=ARM_B_DIR/tumor.markdup.bam
       normal=ARM_B_DIR/normal.markdup.bam
Arm C: Mutect2 tumor-only on ARM_C_DIR/tumor.consensus.bam
       LoFreq tumor-only on same BAM
Arm D: Mutect2 tumor-normal
       tumor=ARM_D_DIR/tumor.consensus.bam
       normal=ARM_D_DIR/normal.consensus.bam

All arms: FilterMutectCalls after Mutect2
Use gnomAD and PON if available, skip gracefully with warning if not
Output: filtered VCFs in respective ARM_X_DIR

### 17_evaluate_results.sh
Core evaluation — most important script for portfolio.

For each arm:
  TRUE POSITIVES:
  Intersect called variants with KNOWN_POSITIVES_VCF
  Stratify by VAF tier (Sample A VAF from INFO field)
  Per tier: TP, FN, sensitivity = TP/(TP+FN)

  FALSE POSITIVES:
  Intersect called variants with KNOWN_NEGATIVES_BED
  FP rate per Mb = FP count / 0.053185

  Per-arm summary: TP, FP, FN, sensitivity, FP rate/Mb, precision, F1

Cross-arm comparisons:
  Arm C vs A: consensus benefit tumor-only (sensitivity gain, FP reduction)
  Arm D vs B: consensus benefit tumor-normal
  Arm B vs A: paired-normal over-suppression (how many variants lost)
  Arm D vs C: value of paired normal with consensus

Output:
  results/evaluation/per_arm_summary.tsv
  results/evaluation/per_variant_calls.tsv
    columns: variant, gene, sample_a_vaf, expected_ef_vaf, tier,
             called_arm_a, called_arm_b, called_arm_c, called_arm_d,
             vaf_arm_a, vaf_arm_b, vaf_arm_c, vaf_arm_d
  results/evaluation/cross_arm_comparison.tsv
  results/evaluation/fp_rate_comparison.tsv

### 18_generate_report.sh
Generate results/REPORT.md covering:
1. Executive summary
2. Dataset and experimental design
3. Methodological note on coordinate consensus vs true UMI calling
4. UMI/consensus metrics (family size distribution, compression ratio)
5. Variant calling results table (all arms, all tiers)
6. Key finding 1: consensus benefit (Arm C vs A)
7. Key finding 2: paired-normal over-suppression (Arm B vs A)
8. Key finding 3: FP rate comparison
9. Limitations and future work (real UMI data, higher VAF tiers)

## Coding Standards

1. Source 00_config.sh as first action after shebang
2. Use require_cmd and require_file from 00_config.sh
3. All scripts idempotent — check output exists, skip with message if so
4. Log steps: echo "[$(date +%Y-%m-%d\ %H:%M:%S)] description"
7. set -euo pipefail (inherited from config)
5. Use $THREADS and $JAVA_XMX throughout
6. Comments explain biological rationale not just commands
7. No hardcoded paths — everything through config variables
8. Temp files to $TMP_DIR, cleaned up on exit

## Instructions for This Session

1. Read ALL files in scripts/ directory first
2. Summarize each existing script — what it does, any issues
3. Propose exact changes to existing scripts
4. Show each new script before writing to disk
5. Wait for explicit approval before writing any file
6. Proceed one script at a time
7. After all scripts exist, propose run order and parallelization opportunities

DO NOT write any code until existing scripts are read and plan is approved.
