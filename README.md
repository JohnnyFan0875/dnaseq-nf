# dnaseq-nf

## Install

1. Install Nextflow (>=24.10.5) from [official website](https://www.nextflow.io/docs/latest/install.html#install-nextflow)
2. Install Docker (>=28.0.1) from [official website](https://docs.docker.com/desktop/)
3. Install required package
   ```bash
   sudo apt install bwa
   ```
4. Create Docker Image
5. Create Reference Files

   1. GATK Resource Bundle files

      ```bash
      bash scripts/download_gatk_resource_bundle.sh
      ```

      Files include:

      1. Homo_sapiens_assembly38.fasta
      2. Homo_sapiens_assembly38.fasta.fai
      3. Homo_sapiens_assembly38.dict
      4. 1000G_omni2.5.hg38.vcf.gz
      5. 1000G_omni2.5.hg38.vcf.gz.tbi
      6. 1000G_phase1.snps.high_confidence.hg38.vcf.gz
      7. 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
      8. hapmap_3.3.hg38.vcf.gz
      9. hapmap_3.3.hg38.vcf.gz.tbi
      10. Homo_sapiens_assembly38.dbsnp138.vcf
      11. Homo_sapiens_assembly38.dbsnp138.vcf.idx
      12. Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
      13. Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

   2. Reference Index

      ```bash
      bwa index reference/fasta/Homo_sapiens_assembly38.fasta
      ```

   3. Ensembl VEP cache
      ```bash
      bash scripts/download_ensembl_vep_cache.sh
      ```

   4. Clinvar VCF file
      ```bash
      bash scripts/download_clinvar_vcf.sh
      ```

## Usage

nextflow run main.nf --project_name <project_name> -resume

## Reference

### Known Sites for BQSR

- [GATK community (BQSR)](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092-Known-Sites-for-BQSR)
- [GATK Best Practice (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator)
- [GATK Resource Bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)
