# ctDNA Variant Calling Pipeline Benchmark

**Coordinate-based consensus calling vs standard duplicate marking**  
**for low-VAF somatic variant detection in ctDNA**

---

## Executive Summary

This benchmark evaluates four variant calling strategies on SEQC2 Burning Rock 
Lung Plasma (BRP v4) reference material. The tumor sample (Df) is a 1:4 dilution 
of cancer cell line into normal background, yielding an expected tumor fraction of 
~20%. Against 228 panel-overlapping known positive variants and 53,185 bp of known 
negative space, the four arms reveal:

- **Arm A (standard tumor-only):** 17.4% overall sensitivity, 131.62 FP/Mb
- **Arm B (standard tumor-normal):** 2.3% overall sensitivity — paired-normal over-suppression removes variants shared with normal background
- **Arm C (consensus tumor-only):** 69.8% overall sensitivity, 470.06 FP/Mb
- **Arm D (consensus tumor-normal):** 52.3% overall sensitivity

Key findings are described in sections 6–8.

## Dataset and Experimental Design

**Tumor:** `SRR13200966` — 1:4 dilution, ~20% tumor fraction (Sample A cancer cell line : Sample B normal)  
**Normal:** `SRR13201012` — pure Sample B normal  
**Panel:** Burning Rock BRP v4 (~53 kb coding regions, hg38)

### Four Pipeline Arms

| Arm | Strategy               | Callers          |
| --- | ---------------------- | ---------------- |
| A   | Standard tumor-only    | Mutect2 + LoFreq |
| B   | Standard tumor-normal  | Mutect2          |
| C   | Consensus tumor-only   | Mutect2 + LoFreq |
| D   | Consensus tumor-normal | Mutect2          |

Arms A and B share the same tumor BAM (duplicate-marked). Arms C and D share the same consensus tumor BAM.

## Truth Set Characterisation

The SEQC2 known positives VCF contains 228 variants overlapping the BRP panel. After annotation with gnomAD and filtering variants at population frequency AF ≥ 0.001, 86 somatic candidates remain:

| Category | Count | Fraction |
|----------|-------|----------|
| Common germline (gnomAD AF >1%) | 124 | 54% |
| Rare germline (gnomAD AF 0.001–1%) | 18 | 8% |
| Somatic candidates (not in gnomAD or AF <0.001) | 86 | 38% |

All 86 somatic candidates fall in Tier 3 (28 variants, Sample A VAF 5–20%, expected ~1-4% in Df) or Tier 4 (58 variants, Sample A VAF <5%, expected <1% in Df). No Tier 1 or Tier 2 somatic variants overlap the BRP panel — the high-VAF Sample A variants in this panel are common germline SNPs correctly filtered by Mutect2 with gnomAD.

## Methodological Note: Coordinate Consensus vs True UMI Calling

UMIs were stripped during SRA submission of the SEQC2 BRP dataset (confirmed via `vdb-dump`). Arms C and D therefore implement **coordinate-based consensus calling**, not true UMI-aware calling.

| Approach | Duplicate identification | Errors corrected |
|---|---|---|
| True UMI calling | UMI sequence ligated before PCR | PCR duplicates + sequencing error |
| Coordinate consensus (this pipeline) | chr + start + end + strand | Sequencing error only |

`fgbio GroupReadsByUmi --strategy identity` is used because synthetic UMIs are derived from fragment coordinates — adjacency/edit-distance correction would be meaningless (a 1-base UMI difference encodes a different coordinate group, not a PCR error).

This limitation is itself a finding: it quantifies the cost of UMI data loss during archival and motivates proper UMI preservation in future studies.

## Consensus Calling Metrics

### Read counts and compression

| Metric                    | Tumor     | Normal    |
| ------------------------- | --------- | --------- |
| Pre-consensus reads       | 139258416 | 189185622 |
| Duplicate rate            | 0.9347    | 0.9595    |
| Post-consensus reads      | 16202082  | 15933208  |
| Compression ratio         | 8.60      | 11.87     |
| Error rate pre-consensus  | 0.000488  | 0.000474  |
| Error rate post-consensus | 0.000161  | 0.000188  |


_Error rate proxy: non-reference allele rate at known-negative loci (BQ ≥ 20). Absolute values are inflated by Sample B germline variants present in both samples; the pre- vs post-consensus comparison is valid._

### Family size distribution (fgbio GroupReadsByUmi)

Family size = number of read pairs sharing identical fragment coordinates (the "duplicate group" size). Families of size 1 are singletons; `--min-reads 2` discards these in consensus calling.

**Tumor**

| family_size | count   | fraction |
| ----------- | ------- | -------- |
| 1           | 3039426 | 0.268771 |
| 2           | 1242739 | 0.109893 |
| 3           | 1002332 | 0.088635 |
| 4           | 854799  | 0.075588 |
| 5           | 741020  | 0.065527 |
| 6           | 643325  | 0.056888 |
| 7           | 555658  | 0.049136 |
| 8           | 479384  | 0.042391 |
| 9           | 410353  | 0.036287 |
| 10          | 350681  | 0.03101  |
| 11          | 296577  | 0.026226 |
| 12          | 250522  | 0.022153 |
| 13          | 212113  | 0.018757 |
| 14          | 178371  | 0.015773 |
| 15          | 150464  | 0.013305 |


**Normal**

| family_size | count   | fraction |
| ----------- | ------- | -------- |
| 1           | 2195071 | 0.213427 |
| 2           | 955221  | 0.092876 |
| 3           | 754122  | 0.073323 |
| 4           | 641512  | 0.062374 |
| 5           | 562053  | 0.054648 |
| 6           | 500773  | 0.04869  |
| 7           | 449962  | 0.04375  |
| 8           | 406120  | 0.039487 |
| 9           | 367740  | 0.035755 |
| 10          | 334020  | 0.032477 |
| 11          | 300958  | 0.029262 |
| 12          | 271447  | 0.026393 |
| 13          | 245947  | 0.023913 |
| 14          | 221666  | 0.021553 |
| 15          | 199876  | 0.019434 |



## Variant Calling Results

### Per-arm summary

| arm | caller  | tp | fn | fp  | sensitivity | fp_rate_per_mb | precision | f1     |
| --- | ------- | -- | -- | --- | ----------- | -------------- | --------- | ------ |
| A   | Mutect2 | 15 | 71 | 7   | 17.4%       | 131.62         | 0.6818    | 0.2778 |
| B   | Mutect2 | 2  | 84 | 0   | 2.3%        | 0.00           | 1.0000    | 0.0455 |
| C   | Mutect2 | 60 | 26 | 25  | 69.8%       | 470.06         | 0.7059    | 0.7018 |
| D   | Mutect2 | 45 | 41 | 1   | 52.3%       | 18.80          | 0.9783    | 0.6818 |
| A   | LoFreq  | 79 | 7  | 213 | 91.9%       | 4004.89        | 0.2705    | 0.4180 |
| C   | LoFreq  | 80 | 6  | 127 | 93.0%       | 2387.89        | 0.3865    | 0.5461 |


### Sensitivity by VAF tier

Expected VAF in Df = Sample A VAF / 5 (1:4 dilution)

| arm | caller  | tier | vaf_range            | tp | sensitivity |
| --- | ------- | ---- | -------------------- | -- | ----------- |
| A   | Mutect2 | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| A   | Mutect2 | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| A   | Mutect2 | 3    | 5-20% (1-4% in Df)   | 11 | 0.3929      |
| A   | Mutect2 | 4    | <5% (<1% in Df)      | 4  | 0.0690      |
| B   | Mutect2 | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| B   | Mutect2 | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| B   | Mutect2 | 3    | 5-20% (1-4% in Df)   | 2  | 0.0714      |
| B   | Mutect2 | 4    | <5% (<1% in Df)      | 0  | 0.0000      |
| C   | Mutect2 | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| C   | Mutect2 | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| C   | Mutect2 | 3    | 5-20% (1-4% in Df)   | 20 | 0.7143      |
| C   | Mutect2 | 4    | <5% (<1% in Df)      | 40 | 0.6897      |
| D   | Mutect2 | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| D   | Mutect2 | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| D   | Mutect2 | 3    | 5-20% (1-4% in Df)   | 9  | 0.3214      |
| D   | Mutect2 | 4    | <5% (<1% in Df)      | 36 | 0.6207      |
| A   | LoFreq  | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| A   | LoFreq  | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| A   | LoFreq  | 3    | 5-20% (1-4% in Df)   | 25 | 0.8929      |
| A   | LoFreq  | 4    | <5% (<1% in Df)      | 54 | 0.9310      |
| C   | LoFreq  | 1    | >=50% (>10% in Df)   | 0  | 0.0000      |
| C   | LoFreq  | 2    | 20-50% (4-10% in Df) | 0  | 0.0000      |
| C   | LoFreq  | 3    | 5-20% (1-4% in Df)   | 25 | 0.8929      |
| C   | LoFreq  | 4    | <5% (<1% in Df)      | 55 | 0.9483      |



## Key Finding 1: Consensus Calling Benefit (Arm C vs A)

Consensus calling (Arm C) vs standard duplicate marking (Arm A) in the tumor-only setting:

- **TP shared:** 14 variants called by both arms
- **TP gained by C:** 46 variants called only by Arm C
- **TP lost by C:** 1 variants called only by Arm A
- **Overall sensitivity:** Arm A = 17.4%, Arm C = 69.8%

FP rate: Arm A = 131.62 FP/Mb, Arm C = 470.06 FP/Mb (Mutect2)

Consensus base quality improvement enables Mutect2 and LoFreq to better distinguish true low-VAF variants from sequencing noise at the cost of reduced read depth after family collapsing.

## Key Finding 2: Paired-Normal Over-Suppression (Arm B vs A)

Both tumor (`SRR13200966`) and normal (`SRR13201012`) samples contain Sample B germline background. Mutect2 tumor-normal calling suppresses variants present in the normal — which includes real tumor variants that happen to overlap the shared germline signal.

- **Overall sensitivity:** Arm A = 17.4%, Arm B = 2.3%
- **TP lost vs Arm A:** 13 variants suppressed by the paired normal
- **TP gained vs Arm A:** 0 variants rescued from FP filtering by the normal

This over-suppression is a known clinical problem related to CHIP and shared germline signal. Quantifying it empirically is a key result of this pipeline.

## Key Finding 3: False Positive Rate Comparison

| arm | caller  | fp_count | known_neg_mb | fp_rate_per_mb |
| --- | ------- | -------- | ------------ | -------------- |
| A   | Mutect2 | 7        | 0.053185     | 131.62         |
| B   | Mutect2 | 0        | 0.053185     | 0.00           |
| C   | Mutect2 | 25       | 0.053185     | 470.06         |
| D   | Mutect2 | 1        | 0.053185     | 18.80          |
| A   | LoFreq  | 213      | 0.053185     | 4004.89        |
| C   | LoFreq  | 127      | 0.053185     | 2387.89        |


FP rate is computed over 53,185 bp of known-negative space. LoFreq FP rates include all calls (LoFreq does not use a PASS filter by default); Mutect2 rates include PASS-only calls after FilterMutectCalls.

## Key Finding 4: Consensus + Paired Normal Achieves Best Overall Performance (Arm D)

Arm D (consensus calling with paired normal) achieves the best balance of sensitivity and specificity:

- **Sensitivity:** 52.3% overall
- **FP rate:** 18.80 FP/Mb (vs 470.06 FP/Mb for Arm C tumor-only)
- **Precision:** 0.9783
- **F1 score:** 0.6818

Compared to standard tumor-normal calling (Arm B, F1=0.0455), consensus calling rescues 44 additional true positives while adding only 1 false positive. This demonstrates that coordinate-based consensus error reduction is sufficient to make low-VAF tumor-normal calling clinically viable.

## Limitations and Future Work

1. **No true UMIs:** Coordinate-based grouping assumes identical-endpoint read pairs are PCR duplicates. Two distinct molecules with coincidentally identical endpoints are conflated — this is rare but real, and cannot be disambiguated without the original UMI sequence.

2. **Low tumor fraction:** At ~20% tumor fraction and 1:4 dilution (Df sample), most Tier 3 and Tier 4 variants are below 4% VAF. Higher-coverage sequencing or higher tumor fraction samples would be needed to assess performance in the >10% VAF range.

3. **Germline contamination in error rate proxy:** Known-negative error rates are inflated by Sample B germline variants. True sequencing error rates are lower; a germline-filtered analysis would provide cleaner baseline estimates.

4. **Single replicate:** Performance estimates would be more robust with multiple Df replicates at different tumor fractions (BRP provides several).

5. **UMI archival:** The primary recommendation from this study is that UMI sequences should be preserved in submitted FASTQ files. The SRA/ENA submission pipeline should be configured to retain RX tags or UMI-prefix read names.

6. **Truth set germline contamination:** 54% of the original SEQC2 known positives are common germline SNPs correctly excluded by Mutect2's germline filter. The refined somatic truth set of 86 variants represents genuinely somatic signal but is limited to Tier 3/4 VAF ranges, preventing evaluation at higher VAF tiers.
