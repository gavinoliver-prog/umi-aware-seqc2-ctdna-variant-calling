# ctDNA Variant Calling Pipeline

## Quick start

```bash
conda env create -f environment.yml
conda activate ctdna2
bash scripts/run_all.sh
```

The pipeline benchmarks coordinate-based consensus calling (fgbio) against
standard duplicate marking (Picard) for low-VAF somatic variant detection in
ctDNA, using SEQC2 BRP reference samples with published ground truth.

See `CLAUDE.md` for full experimental design and methodology.

---

## Note on choice of dataset

This project builds directly on a prior SEQC2-based ctDNA variant-calling workflow in which we evaluated ultra-deep targeted sequencing using the Illumina arm of the SEQC2 liquid biopsy study. That initial work intentionally focused on UMI-unaware analysis and demonstrated expected limitations in low-VAF detection despite increasing depth. In extending this work to a UMI-aware framework, we attempted to reuse SEQC2 data across multiple assay arms (including Burning Rock and Roche) under the assumption that molecular barcodes would be accessible in the public data. However, inspection of both SRA-derived FASTQs and BAMs showed that UMIs were not recoverable in a usable form (e.g., absent from headers or missing from standard tags such as RX), despite being described in the original study. Given this limitation, we pivoted to a QIAseq/smCounter-style dataset with verifiable, accessible UMIs, allowing us to implement and rigorously evaluate UMI-aware processing while preserving the overall philosophy of using real data, real tools, and reproducible workflows.
