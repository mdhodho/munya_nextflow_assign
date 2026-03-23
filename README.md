#Reproducible Bioinformatics Pipeline (Nextflow + Singularity)

## Overview

This project implements a reproducible Nextflow workflow for processing raw sequencing data from FASTQ files to variant calling (VCF). The pipeline integrates containerization (Singularity) and workflow automation (Nextflow).

## Workflow

FASTQ → FastQC → Trimmomatic → FastQC → BWA → BAM → Variant Calling → VCF/

## Data Source

Raw sequencing data were obtained from the European Nucleotide Archive (ENA):

* Run accession: **ERR1955323**

Download:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_2.fastq.gz
```

Rename:

```bash
mv ERR1955323_1.fastq.gz sample1_R1.fastq.gz
mv ERR1955323_2.fastq.gz sample1_R2.fastq.gz
```

## Tools Used

* FastQC
* Trimmomatic (Singularity container)
* BWA-MEM
* Samtools
* Bcftools
* Nextflow

## Containerization

```bash
singularity build --fakeroot trimmomatic.sif trimmomatic.def
```


## ▶️ Running the Pipeline

```bash
./nextflow run main.nf --threads 4
```

Resume:

```bash
./nextflow run main.nf --threads 4 -resume
```

## Inputs

* data/sample1_R1.fastq.gz
* data/sample1_R2.fastq.gz
* ref/reference.fa
* ref/adapters.fa

## 📁 Outputs

* FastQC reports
* Trimmed reads
* BAM files
* Final VCF:

results/vcf/sample1.vcf

## 📊 Results

A total of **2,724 variants** were identified in the sample.

## Key Concepts

### Reference Genome

Used to align reads and detect variants.

### Adapter Sequences

Removed to improve alignment accuracy.

## Reproducibility

* Nextflow workflow
* Singularity container
* Resume capability (-resume)


## Author

Munyaradzi Dhodho

