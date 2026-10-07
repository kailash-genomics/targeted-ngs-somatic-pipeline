# Targeted NGS Somatic Variant Calling Pipeline

An end-to-end **Snakemake-based pipeline for tumor-only targeted NGS somatic variant analysis using hg19/GRCh37**.

The pipeline processes paired-end FASTQ files through quality control, read preprocessing, alignment, BAM processing, target-region coverage analysis, somatic variant calling, variant filtering, annotation, and automated QC reporting.

**Scope:** SNV/indel calling, variant filtering, annotation, target coverage metrics, and automated QC worksheet generation.

> **Not implemented:** CNV, fusion, MSI, TMB, matched tumor-normal analysis, and germline variant analysis.

---

## Workflow

```text
Paired-end FASTQ
       |
       v
  FastQC / fastp
       |
       v
    BWA-MEM
       |
       v
 BAM Processing
       |
       +----------------------+
       |                      |
       v                      v
     Sort              MarkDuplicates
       |                      |
       +----------+-----------+
                  |
                  v
         Target Coverage QC
                  |
                  v
             GATK Mutect2
                  |
       +----------+-----------+
       |          |           |
       v          v           v
 LearnRead   GetPileup    Calculate
 Orientation  Summaries   Contamination
 Model
       |          |           |
       +----------+-----------+
                  |
                  v
         FilterMutectCalls
                  |
                  v
          VCF Normalization
                  |
                  v
     SnpEff / ClinVar / CIViC
                  |
                  v
      Variant + Coverage Tables
                  |
                  v
       Automated QC Worksheet
```

---

## Pipeline Components

| Stage                         | Tools                                            |
| ----------------------------- | ------------------------------------------------ |
| FASTQ quality control         | FastQC                                           |
| Read preprocessing            | fastp                                            |
| Alignment                     | BWA-MEM                                          |
| BAM processing                | SAMtools / GATK                                  |
| Duplicate marking             | GATK MarkDuplicates                              |
| Coverage analysis             | mosdepth / bedtools / GATK                       |
| Somatic variant calling       | GATK Mutect2                                     |
| Orientation-bias modeling     | GATK LearnReadOrientationModel                   |
| Contamination estimation      | GATK GetPileupSummaries / CalculateContamination |
| Variant filtering             | GATK FilterMutectCalls                           |
| VCF processing                | BCFtools / bgzip / tabix                         |
| Variant annotation            | SnpEff / ClinVar                                 |
| Clinical evidence integration | CIViC                                            |
| Transcript prioritization     | MANE                                             |
| QC and reporting              | Python / python-docx                             |
| Workflow management           | Snakemake                                        |

---

## Analysis Workflow

### 1. FASTQ Quality Control

Raw paired-end FASTQ files are processed using **fastp** for adapter and quality trimming.

**FastQC** is used for FASTQ quality assessment.

The current workflow runs FastQC on the R1 input and performs paired-end preprocessing with fastp.

### 2. Alignment

Trimmed paired-end reads are aligned to the **hg19/GRCh37 reference genome** using **BWA-MEM**.

Read-group information is added during alignment.

### 3. BAM Processing

The pipeline performs:

* BAM sorting
* BAM indexing
* Duplicate marking using GATK MarkDuplicates
* Alignment statistics generation
* Target-region hybrid-selection metrics

### 4. Target Coverage Analysis

Target-region coverage is calculated using the configured panel BED file.

The workflow generates coverage metrics at multiple depth thresholds, including:

* ≥100X
* ≥250X
* ≥500X
* ≥1000X

Gene-level and target-region coverage tables are generated for QC assessment.

### 5. Tumor-Only Somatic Variant Calling

Somatic SNVs and indels are called using **GATK Mutect2** in a tumor-only workflow.

The workflow additionally performs:

* Read-orientation bias modeling using `LearnReadOrientationModel`
* Pileup-based estimation using `GetPileupSummaries`
* Contamination estimation using `CalculateContamination`
* Final variant filtering using `FilterMutectCalls`

A germline population allele-frequency resource is used as part of the Mutect2 workflow.

### 6. VCF Normalization and Processing

Variant calls undergo normalization and processing using **BCFtools**.

Compressed VCF files are indexed using **tabix**.

### 7. Variant Annotation

Variants are annotated using **SnpEff**.

The reporting workflow integrates additional annotation and evidence resources including:

* **ClinVar**
* **CIViC**
* **MANE transcript information**
* Population allele-frequency information where available

The final variant table contains annotation, sequencing-support, and review-related fields for downstream assessment.

### 8. Coverage and Variant Reporting

The workflow generates tables containing:

* Variant-level information
* Read depth
* Alternate read count
* Variant allele fraction (VAF)
* Annotation information
* Population allele-frequency information
* ClinVar information
* Coverage metrics

### 9. Automated QC Reporting

The pipeline generates a structured **Word QC worksheet** containing sequencing, alignment, coverage, variant, and QC metrics.

A provenance JSON file is also generated to record relevant input, database, and pipeline information used during analysis.

---

## Requirements

### Operating System

The workflow is designed for:

* Linux
* WSL2
* macOS

### Conda Environment

The main workflow is designed to run in a Conda environment named:

```bash
somatic
```

Required command-line tools include:

```text
fastp
fastqc
bwa
samtools
bcftools
tabix
bgzip
gatk
mosdepth
bedtools
snakemake
java
```

SnpEff is configured through a separate environment:

```text
snpeff_env
```

Required Python packages include:

```text
pandas
pyyaml
python-docx
pysam
```

> Exact software and database versions should be recorded when reproducing an analysis.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/kailash-genomics/targeted-ngs-somatic-pipeline.git
cd targeted-ngs-somatic-pipeline
```

### 2. Activate the Conda environment

```bash
conda activate somatic
```

### 3. Configure reference files and databases

The large reference files and databases are intentionally excluded from this public repository.

Place or symlink the required resources locally:

```text
ref/
├── hg19.fa
└── hg19.fa.*              # BWA/reference indexes

bed/
└── target_regions.hg19.bed

db/
├── clinvar/
├── mutect2/
├── civic/
└── other required resources
```

The exact resource paths are configured in:

```text
config/settings.yaml
```

### 4. Configure the pipeline

Edit:

```text
config/settings.yaml
```

Configure:

* Reference genome
* Target BED file
* Annotation/database paths
* Number of threads
* Java memory
* QC thresholds
* Variant filtering thresholds
* Panel name
* Other analysis parameters

The target BED path is configurable, allowing the workflow to be adapted to different targeted panels.

### 5. Run preflight checks

Run the standard preflight checks:

```bash
./run_sample.sh check
```

For additional checks:

```bash
./run_sample.sh check --full
```

The preflight workflow checks important inputs, tools, reference indexes, databases, BED-file properties, computational resources, and other pipeline requirements.

### 6. Run a sample

```bash
./run_sample.sh SAMPLE_ID /path/to/R1.fastq.gz /path/to/R2.fastq.gz
```

### 7. Check sample status

```bash
./run_sample.sh status SAMPLE_ID
```

---

## Repository Structure

```text
targeted-ngs-somatic-pipeline/
│
├── Snakefile
├── run_sample.sh
│
├── config/
│   ├── metadata.csv
│   ├── samples.csv
│   ├── settings.yaml
│   └── template.docx
│
├── scripts/
│   ├── bamstats.sh
│   ├── coverage_tables.py
│   ├── fill_report.py
│   ├── final_table.py
│   └── preflight.py
│
├── .gitignore
├── LICENSE
└── README.md
```

Large sequencing files, reference genomes, databases, intermediate results, and runtime files are excluded from the public repository through `.gitignore`.

---

## Output

For a completed sample, the final QC worksheet is copied to:

```text
final_reports/SAMPLE_ID/SAMPLE_ID.QC_worksheet.docx
```

The workflow also generates intermediate files required for:

* FASTQ preprocessing
* Alignment
* BAM processing
* QC
* Coverage analysis
* Somatic variant calling
* Variant filtering
* Annotation
* Variant tables

A provenance JSON file is generated alongside the final reporting outputs.

---

## Configuration

Pipeline parameters are managed through:

```text
config/settings.yaml
```

The configuration includes:

* Reference genome
* Target BED file
* Database/resource locations
* Computational resources
* QC thresholds
* Variant support thresholds
* Panel name
* Annotation resources

Example thresholds are provided for development and demonstration purposes.

**Important:** QC and variant-filtering thresholds must be evaluated and analytically validated for the intended assay, panel, sequencing platform, and laboratory workflow before clinical use.

---

## Data Privacy

Do **not** commit patient-identifiable or confidential clinical data to this public repository.

Do not upload:

* Patient FASTQ files
* BAM/CRAM files
* Patient VCF files
* Clinical reports
* HPE reports
* Patient identifiers
* Hospital-confidential information
* Passwords or API credentials
* Internal laboratory data

Use synthetic, publicly available, or appropriately de-identified data for demonstrations.

---

## Clinical Disclaimer

This repository contains **research/development software** and is not a clinically validated diagnostic pipeline.

The included QC and variant-filtering thresholds are example values and require appropriate analytical validation before use in a clinical setting.

**Human review is required before any clinical interpretation.**

Clinical implementation requires appropriate laboratory validation, quality assurance, performance assessment, documentation, and applicable regulatory/accreditation requirements.

---

## Current Scope

### Implemented

* [x] FASTQ quality control
* [x] Read preprocessing
* [x] BWA-MEM alignment
* [x] BAM sorting and indexing
* [x] Duplicate marking
* [x] Alignment statistics
* [x] Target-region coverage analysis
* [x] Tumor-only Mutect2 variant calling
* [x] Read-orientation bias modeling
* [x] Contamination estimation
* [x] Mutect2 variant filtering
* [x] VCF normalization
* [x] SnpEff annotation
* [x] ClinVar integration
* [x] CIViC evidence integration
* [x] MANE transcript information
* [x] Population allele-frequency information
* [x] Variant-level reporting tables
* [x] Coverage QC tables
* [x] Automated Word QC worksheet
* [x] Provenance information generation
* [x] Pipeline preflight validation

### Not Implemented

* [ ] CNV analysis
* [ ] Fusion detection
* [ ] MSI analysis
* [ ] TMB analysis
* [ ] Matched tumor-normal analysis
* [ ] Germline variant analysis

---

## Future Development

Potential future improvements include:

* Expanded variant annotation resources
* COSMIC integration
* Improved coverage visualization
* MultiQC integration
* CNV analysis
* MSI analysis
* TMB estimation
* Fusion detection
* Matched tumor-normal analysis
* Containerized execution
* Automated variant prioritization
* Expanded report generation
* Reproducible test/demo dataset
* Automated workflow testing
* Environment/version locking

---

## Limitations

This workflow is designed as a **research/development pipeline** and has several limitations:

* It currently performs **tumor-only** somatic variant calling.
* Matched tumor-normal analysis is not implemented.
* CNV, fusion, MSI, and TMB analysis are not implemented.
* The configured QC and variant-filtering thresholds are example values and are not clinical validation thresholds.
* Final variant interpretation requires appropriate expert review.
* Reference genomes, annotation databases, population resources, and target BED files must be supplied separately.
* Results depend on the reference genome, panel design, sequencing platform, database versions, and analysis configuration.

---

## Reproducibility and Provenance

The pipeline is designed to retain analysis provenance through generated metadata and provenance information.

Relevant analysis inputs include:

* Sample information
* FASTQ files
* Reference genome
* Target BED file
* Annotation resources
* Population resources
* ClinVar/CIViC resources
* Pipeline configuration
* Tool versions
* Database information

For reproducible analysis, users should record the exact versions and source dates of reference genomes, annotation databases, and software used.

---

## Author

**Kailash Saini**

NGS & Clinical Genomics | Molecular Diagnostics | Bioinformatics

GitHub: [kailash-genomics](https://github.com/kailash-genomics)

---

## License

See the [LICENSE](LICENSE) file for licensing information.
