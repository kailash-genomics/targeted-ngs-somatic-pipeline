# Targeted NGS Somatic Variant Calling Pipeline

An end-to-end **Snakemake-based pipeline for tumor-only targeted NGS somatic variant analysis using hg19/GRCh37**.

The pipeline processes targeted sequencing data from **paired-end FASTQ files** through quality control, read preprocessing, alignment, BAM processing, target coverage analysis, somatic variant calling, filtering, annotation, and automated QC reporting.

**Scope:** SNV/indel calling, variant filtering, annotation, target coverage metrics, and automated QC worksheet generation.

> **Not implemented:** CNV, fusion, MSI, TMB, matched tumor-normal analysis, and germline variant analysis.

---

## 🧬 Workflow

```text
Paired-end FASTQ
       │
       ▼
   FastQC / fastp
       │
       ▼
     BWA-MEM
       │
       ▼
 BAM Processing
 ┌─────┼─────────┐
 │     │         │
Sort  Index   MarkDuplicates
 └─────┼─────────┘
       │
       ▼
 Target Coverage QC
       │
       ▼
    GATK Mutect2
       │
       ▼
 Variant Filtering
       │
       ▼
 VCF Normalization
       │
       ▼
 SnpEff / SnpSift
       │
       ▼
 QC + Final Report
```

---

## 🔬 Pipeline Components

| Stage                   | Tools                |
| ----------------------- | -------------------- |
| FASTQ QC                | FastQC               |
| Read preprocessing      | fastp                |
| Alignment               | BWA-MEM              |
| BAM processing          | SAMtools / GATK      |
| Duplicate marking       | GATK                 |
| Coverage analysis       | mosdepth / bedtools  |
| Somatic variant calling | GATK Mutect2         |
| Variant processing      | BCFtools             |
| Variant annotation      | SnpEff / SnpSift     |
| QC & reporting          | Python / python-docx |
| Workflow management     | Snakemake            |

---

## 📊 Analysis Workflow

### 1. FASTQ Quality Control

Raw paired-end FASTQ files are evaluated using **FastQC** and processed using **fastp** for adapter and quality trimming.

### 2. Alignment

Trimmed reads are aligned to the **hg19/GRCh37 reference genome** using **BWA-MEM**.

### 3. BAM Processing

The pipeline performs:

* BAM sorting
* BAM indexing
* Duplicate marking
* Alignment statistics generation

### 4. Target Coverage Analysis

Target-region coverage is calculated using the configured panel BED file.

The workflow generates coverage metrics used for sequencing/QC assessment, including depth-based metrics.

### 5. Somatic Variant Calling

Somatic SNVs and indels are called using **GATK Mutect2** in a tumor-only workflow.

### 6. Variant Filtering & Normalization

Variant calls undergo configured filtering and VCF normalization using **BCFtools** and related tools.

### 7. Variant Annotation

Variants are annotated using **SnpEff/SnpSift** for downstream analysis and reporting.

### 8. Automated QC Reporting

The pipeline generates a structured **Word QC worksheet** containing relevant QC and analysis metrics.

---

## ⚙️ Requirements

### Operating System

* Linux
* WSL2
* macOS

### Conda Environment

The pipeline is designed to run in a Conda environment named:

```bash
somatic
```

Required tools include:

```text
fastp
fastqc
bwa
samtools
bcftools
GATK >= 4.6
mosdepth
bedtools
snakemake
snpEff
```

Required Python packages:

```text
pandas
pyyaml
python-docx
pysam
```

---

## 🚀 Quick Start

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

Place or symlink the required resources:

```text
ref/
├── hg19.fa
└── hg19.fa.*        # reference indexes

bed/
└── target_regions.hg19.bed

db/
└── required databases/resources
```

### 4. Configure the pipeline

Edit:

```text
config/settings.yaml
```

Configure:

* Reference genome
* Target BED file
* Database/resource paths
* Number of threads
* Memory
* QC thresholds
* Other analysis parameters

The BED file path is configurable, allowing the workflow to be adapted to different targeted panels.

### 5. Run preflight checks

```bash
./run_sample.sh check
```

### 6. Run a sample

```bash
./run_sample.sh SAMPLE_ID /path/to/R1.fastq.gz /path/to/R2.fastq.gz
```

### 7. Check sample status

```bash
./run_sample.sh status SAMPLE_ID
```

---

## 📁 Repository Structure

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

---

## 📄 Output

The final QC worksheet is generated at:

```text
final_reports/SAMPLE_ID/SAMPLE_ID.QC_worksheet.docx
```

The pipeline also generates intermediate files required for alignment, QC, coverage analysis, variant calling, and annotation.

---

## ⚙️ Configuration

Pipeline parameters are managed through:

```text
config/settings.yaml
```

This includes configurable settings for:

* Reference genome
* Target BED file
* Database/resource locations
* Computational resources
* QC thresholds
* Pipeline parameters

**Important:** QC and variant filtering thresholds provided in the configuration are example values and must be evaluated and validated for the intended assay and sequencing workflow.

---

## 🔒 Data Privacy

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

Use synthetic or appropriately de-identified data for demonstrations.

---

## ⚠️ Clinical Disclaimer

This repository contains **research/development software** and is not a clinically validated diagnostic pipeline.

The included QC and variant-filtering thresholds are examples and require appropriate analytical validation before use in a clinical setting.

**Human review is required before any clinical interpretation.**

Clinical implementation requires appropriate laboratory validation, quality assurance, performance assessment, documentation, and applicable regulatory/accreditation requirements.

---

## 🚧 Current Scope

### Implemented

* [x] FASTQ quality control
* [x] Read preprocessing
* [x] BWA-MEM alignment
* [x] BAM processing
* [x] Duplicate marking
* [x] Target-region coverage analysis
* [x] Tumor-only Mutect2 variant calling
* [x] Variant filtering
* [x] VCF normalization
* [x] SnpEff/SnpSift annotation
* [x] QC metric generation
* [x] Automated Word QC worksheet

### Not Implemented

* [ ] CNV analysis
* [ ] Fusion detection
* [ ] MSI analysis
* [ ] TMB analysis
* [ ] Matched tumor-normal analysis
* [ ] Germline variant analysis

---

## 🔮 Future Development

Potential future improvements include:

* Expanded variant annotation resources
* ClinVar integration
* COSMIC annotation
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

---

## 👨‍🔬 Author

**Kailash Saini**

NGS & Clinical Genomics | Molecular Diagnostics | Bioinformatics

GitHub: [kailash-genomics](https://github.com/kailash-genomics)

---

## 📜 License

See the [LICENSE](LICENSE) file for licensing information.
