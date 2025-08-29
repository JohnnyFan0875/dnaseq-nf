# dnaseq-nf

## Introduction

**dnaseq-nf** is a Nextflow-based pipeline designed for end-to-end DNA sequencing analysis. It automates alignment, variant calling, and annotation steps following GATK Best Practices. The pipeline integrates widely used bioinformatics tools (BWA, GATK, VEP, ClinVar) and public reference datasets, enabling reproducible and scalable analysis of whole-genome sequencing (WGS) data.

Key features:

- Fully containerized workflow (Docker)
- Implements GATK Best Practices for variant calling
- Automated download of reference resources (GATK bundle, VEP cache, ClinVar, TaiwanGenomes)
- Generates high-quality variant calls annotated with population and clinical databases

## Workflow

The pipeline consists of the following steps:

1. **Preprocessing**

   - Read alignment with BWA-MEM
   - Sorting and marking duplicates
   - Base Quality Score Recalibration (BQSR) using GATK and known variant sites

2. **Variant Calling**

   - SNP and indel discovery with GATK HaplotypeCaller

3. **Annotation**

   - Functional annotation using Ensembl VEP
   - Clinical annotation with ClinVar
   - Population frequency annotation (TaiwanGenomes)

4. **Results**

   - High-confidence VCF files with annotations
   - Quality control metrics for alignment and variant calling

## Installation

1. Install Nextflow (>=24.10.5) from [official website](https://www.nextflow.io/docs/latest/install.html#install-nextflow)
2. Install Docker (>=28.0.1) from [official website](https://docs.docker.com/desktop/)
3. Install required packages:

   ```bash
   sudo apt install bwa rsync tar wget
   ```

4. Download Reference Files

   1. GATK Resource Bundle files

      ```bash
      bash scripts/download_gatk_resource_bundle.sh
      ```

   2. Ensembl VEP cache

      ```bash
      bash scripts/download_ensembl_vep_cache.sh
      ```

   3. Clinvar VCF file

      ```bash
      bash scripts/download_clinvar_vcf.sh
      ```

   4. TaiwanGenome VCF file

      ```bash
      bash scripts/download_taiwangenome.sh
      ```

## Usage

```bash
nextflow run main.nf --project_name <project_name> -resume
```

Example:

```bash
nextflow run main.nf --project_name my_WGS_analysis -resume
```

## References

- [GATK community (Known Sites for BQSR)](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092-Known-Sites-for-BQSR)
- [GATK Best Practice (VariantRecalibrator)](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator)
- [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)
