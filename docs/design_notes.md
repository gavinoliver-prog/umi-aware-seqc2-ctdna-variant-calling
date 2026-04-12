Dataset Selection and Experimental Design: Decision Log
Starting point: Why SEQC2?
The SEQC2 (Sequencing Quality Control Phase 2) consortium published a multi-site, cross-platform ctDNA benchmarking study using purpose-built reference samples with known variants at defined allele frequencies. This makes it one of the few publicly available datasets where:

Ground truth variants are independently validated (including by ddPCR)
Multiple tumor fractions are available from the same source material
Multiple sequencing platforms and panel designs were applied to the same samples
Known negative regions are explicitly defined, enabling false positive rate calculation

This combination of features — real biology, known truth, multiple VAF tiers, explicit negative space — is rare in public data and makes SEQC2 uniquely suitable for benchmarking variant calling pipelines.
Why the BRP (Burning Rock) assay arm?
The SEQC2 study included five ctDNA assay platforms. The Burning Rock Lung Plasma v4 (BRP) arm was selected for the following reasons:
Panel size and enrichment efficiency. BRP targets approximately 50kb across ~70 clinically relevant oncology genes using an amplicon-based approach. Amplicon panels achieve higher on-target rates (>80%) than hybrid capture panels, meaning a larger fraction of sequenced reads contribute to variant calling. The Illumina TruSight Tumor 170 (TST170) arm, initially considered, targets ~500kb with hybrid capture and showed ~25% on-target rate in this dataset — yielding insufficient depth per target base for meaningful UMI family formation.
Sequencing depth. The BRP Ef replicate (SRR13200999) contains approximately 58.5 million read pairs. At >80% on-target rate across a 50kb panel this yields an estimated raw depth of ~140,000x, with an expected duplicate rate of ~93% at this depth — corresponding to approximately 10,000x unique molecule coverage and mean UMI family sizes of ~14 reads. This is the regime where UMI consensus calling is demonstrably beneficial.
UMI library preparation. The BRP assay used the Burning Rock HS UMI library preparation kit, meaning UMIs were present in the original experiment. Although UMI sequences were not preserved in the SRA submission (a systematic issue with SRA's ETL normalization pipeline, confirmed via vdb-dump --info showing FMT: FASTQ and LDR: latf-load), the underlying duplicate structure from the original UMI-tagged library preparation is preserved in fragment coordinates.
The UMI availability problem
A significant portion of the project scoping involved investigating whether UMI sequences were accessible in public data. The findings were:
SRA strips UMIs systematically. When sequencing data is submitted to NCBI SRA as FASTQ files (after upstream lab processing), any UMI sequences that were trimmed from reads before submission are permanently lost. The SRA ETL normalization pipeline does not recover them. This was confirmed for multiple datasets:

SEQC2 ILM arm (SRR19193668): vdb-dump showed FMT: FASTQ, LDR: latf-load — UMIs absent
SEQC2 BRP arm (SRR13201012, SRR13200999): ENA filereport showed empty submitted_ftp field with only standardized FASTQs available — original BAM not deposited
PRJNA844028 (structured UMI paper): Read length distribution showed biological insert sizes, not fixed UMI-prefixed reads — UMIs stripped pre-submission

Original submitted files are not always recoverable. NCBI maintains only ETL-normalized data online. Original submitted BAMs (which would contain RX:Z: UMI tags) are stored only in AWS/GCP cloud buckets and are not always present even there. ENA's submitted_ftp field is the most reliable indicator — an empty field means the original format was FASTQ and UMIs are unrecoverable.
The solution: fragment coordinate-based UMI reconstruction. For cfDNA data, fragment endpoints are biologically determined by nucleosomal protection during apoptotic DNA release. Read pairs sharing identical chromosome, start position, end position, and strand orientation originate from the same original DNA molecule — the same relationship that UMIs encode. Assigning deterministic synthetic UMIs based on a hash of fragment coordinates reconstructs the molecule-family structure from the authentic PCR duplicate patterns present in the data, without fabricating any biological signal.
Sample selection: Ef as tumor, Bf as normal
The SEQC2 BRP dataset includes multiple samples representing different tumor fractions:
Bf  (SRR13201012) — 0:1 tumor:normal   — pure normal background
Df  (SRR13200999) — 1:4 tumor:normal   — ~20% tumor fraction (LBx-high)
Ef  (SRR13200999) — 1:24 tumor:normal  — ~4% tumor fraction (LBx-low)
Ff                — 1:124 tumor:normal — ~0.8% tumor fraction
Ef was selected as the tumor sample because the ~4% tumor fraction places known variants at approximately 0.16–4% VAF in the sequenced sample (Sample A VAF ÷ 25 dilution factor). This range spans the clinically relevant ctDNA detection zone for early-stage cancer and MRD monitoring, and represents the regime where UMI consensus calling provides the greatest benefit over standard duplicate marking.
Bf was selected as the paired normal because it represents the pure Sample B background with zero tumor fraction, matching the normal cell line component of the Ef mixture.
Truth set construction
Known positives. The SEQC2 consortium published 42,317 known positive variants for Sample A in hg38 coordinates (lifted over from hg19, with 39 variants failing liftover). These represent variants present in the Sample A cancer cell line mixture at defined VAFs, validated across multiple sequencing platforms.
The BRP panel BED file was originally in hg19 coordinates and required liftover to hg38 using UCSC liftOver before intersection with the known positives VCF. After liftover, intersection of the known positives with the BRP panel identified 228 evaluable variants across 70 genes.
These 228 variants are stratified into four VAF tiers based on their frequency in pure Sample A (expected VAF in Ef = Sample A VAF ÷ 25):
TierSample A VAFExpected Ef VAFVariant countInterpretation1≥50%>2%30Detectable by standard calling220–50%0.8–2%39Borderline; UMI should help35–20%0.2–0.8%76UMI arms expected to outperform4<5%<0.2%83Below practical detection limit
Key genes represented include EGFR, TP53, BRCA1, BRCA2, KRAS, BRAF, ALK, PIK3CG, KEAP1, STK11, CDKN2A, APC, and RET — a clinically comprehensive oncology gene set.
Known negatives. The SEQC2 consortium also published known negative regions — genomic positions confirmed to have no variants in Sample A. Intersection of the known negatives BED with the BRP panel identified 53,185 bp of evaluable negative space. Any variant called within this space is by definition a false positive, enabling calculation of false positive rate expressed as variants per megabase of known negative space — the standard metric used in clinical NGS validation.
Experimental design summary
Four-arm parallel comparison:

Arm A: Standard tumor-only
       Ef → BWA-MEM → Picard MarkDuplicates → Mutect2
       
Arm B: Standard tumor-normal
       Ef+Bf → BWA-MEM → Picard MarkDuplicates → Mutect2
       Expected finding: paired-normal over-suppression due to shared
       germline variants between Ef and Bf (both contain Sample B background)
       
Arm C: UMI-aware tumor-only  
       Ef → BWA-MEM → fgbio UMI reconstruction → Mutect2
       
Arm D: UMI-aware tumor-normal
       Ef+Bf → BWA-MEM → fgbio UMI reconstruction → Mutect2

Primary metrics:
- Sensitivity per VAF tier (TP / known positives per tier)
- False positive rate per Mb of known negative space
- F1 score per arm

Key comparisons:
- Arm C vs A: UMI benefit in tumor-only setting
- Arm D vs B: UMI benefit in tumor-normal setting
- Arm B vs A: paired-normal over-suppression (variants lost due to shared background)
- Arm D vs C: added value of paired normal when using UMI consensus

### Final Sample Selection and Accession Mapping

After initial analysis with Ef (1:24 dilution, ~4% tumor fraction), all evaluable 
somatic variants were found to be at <0.8% expected VAF in the sequenced sample — 
below the reliable detection threshold for standard somatic calling. Df (1:4 dilution, 
~20% tumor fraction) was selected as the tumor sample to provide detectable somatic 
signal across multiple VAF tiers.

Site matching was prioritised: Bf normal (Site 25) is paired with Df Site 25 LIB1 
to minimise site-specific technical variation in tumor-normal calling.

| Sample | Role | Dilution | Tumor Fraction | SAMN | SRX | SRR | Library | Site |
|--------|------|----------|----------------|------|-----|-----|---------|------|
| Df | Tumor | 1:4 | ~20% | SAMN16786374 | SRX9634502 | SRR13200966 | SampleDf_BRP2_ST25_25ng_LIB1 | 25 |
| Bf | Normal | 0:1 | 0% | SAMN16786376 | - | SRR13201012 | SampleBf_BRP2_ST25_25ng_LIB1 | 25 |

Previously evaluated (not used in final pipeline):
| Ef | - | 1:24 | ~4% | SAMN16786375 | SRX9634460 | SRR13200999 | SampleEf_BRP2_ST25_25ng_LIB1 | 25 |
| Df (Site 26) | - | 1:4 | ~20% | SAMN16786374 | SRX9634460 | SRR13201008 | SampleDf_BRP2_ST26_25ng_LIB1 | 26 |

### Expected VAF in Df

Expected VAF in Df = Sample A VAF / 5 (1:4 dilution factor)

| Tier | Sample A VAF | Expected Df VAF | Variant count | Interpretation |
|------|-------------|-----------------|---------------|----------------|
| 1 | ≥50% | ≥10% | 30 | Easily detectable by standard calling |
| 2 | 20–50% | 4–10% | 39 | Detectable by standard calling |
| 3 | 5–20% | 1–4% | 76 | Detectable with good sensitivity |
| 4 | <5% | <1% | 83 | Challenging, UMI consensus may help |

Note: VAF tier counts based on original 228-variant truth set. After gnomAD 
germline filtering, 86 somatic candidates remain (28 Tier 3, 58 Tier 4).
At Df dilution, Tier 3 somatic candidates appear at 0.2-4% VAF — within 
the detectable range for standard calling with sufficient depth.

### Truth Set Refinement: Germline Filtering

Initial analysis revealed that the SEQC2 known positives VCF contains 
a mixture of germline and somatic variants. After annotating with gnomAD 
(af-only-gnomad.hg38.vcf.gz) and filtering variants at population 
frequency AF >= 0.001:

| Category | Count | Fraction | Notes |
|----------|-------|----------|-------|
| Common germline (gnomAD AF >1%) | 124 | 54% | Correctly filtered by Mutect2 germline model |
| Rare germline (gnomAD AF 0.001-1%) | 18 | 8% | Partially filtered |
| Somatic candidates (gnomAD AF <0.001 or absent) | 86 | 38% | Retained for evaluation |
| Total | 228 | 100% | |

The 124 common germline variants were correctly identified and filtered 
by Mutect2 using the gnomAD germline resource — this is expected and 
appropriate behaviour. The SEQC2 truth set was designed for broad 
benchmarking across platforms and includes germline variants as 
"known positives" because they are measurable in Sample A. However 
for somatic variant calling evaluation they should be excluded.

The refined somatic-only truth set (86 variants) is stored at:
  ref/SampleA_ref/KnownPositives_hg38.somatic_only.vcf.gz

All 86 somatic candidates fall in Tier 3 (28 variants, Sample A VAF 
5-20%) or Tier 4 (58 variants, Sample A VAF <5%). At Df dilution 
(1:4), expected VAFs are:
  Tier 3: 1-4%  → detectable by standard calling at sufficient depth
  Tier 4: <1%   → challenging, consensus calling expected to help

This VAF range is clinically relevant for ctDNA MRD monitoring where 
tumor fractions of 1-5% are typical in early relapse detection.
