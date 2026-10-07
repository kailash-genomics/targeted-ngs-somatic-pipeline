# Targeted NGS Somatic Variant Pipeline

Tumor-only targeted-panel somatic variant calling pipeline for hg19.

**Scope:** SNV/indel calling, filtering, annotation, coverage metrics, and a Word QC worksheet.
CNV, fusion, MSI, and TMB analysis are **not** implemented.

## Requirements

- Linux / WSL2 / macOS
- Conda environment named `somatic` containing:
  fastp, fastqc, bwa, samtools, bcftools, GATK >=4.6, mosdepth, bedtools, snakemake, snpEff,
  and Python packages (pandas, pyyaml, python-docx, pysam)

## Quick start

1. Clone the repo
2. Place or symlink reference (ref/hg19.fa + indexes), BED (bed/target_regions.hg19.bed), and databases (db/)
3. conda activate somatic
4. ./run_sample.sh check
5. ./run_sample.sh SAMPLE_ID /path/to/R1.fastq.gz /path/to/R2.fastq.gz
6. ./run_sample.sh status SAMPLE_ID

Report: final_reports/SAMPLE_ID/SAMPLE_ID.QC_worksheet.docx

## Configuration

Edit config/settings.yaml for paths, threads, memory, and QC thresholds.
The BED path is configurable — any panel BED can be used.

Thresholds are examples only and must be validated before clinical use.

## Important

- Research / development software
- Unvalidated placeholder thresholds
- Human review required before any clinical interpretation
