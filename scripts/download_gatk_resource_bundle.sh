#!/bin/bash

mkdir reference/fasta
mkdir reference/broad_hg38

# fasta file for alignment
wget -P reference/fasta/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget -P reference/fasta/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget -P reference/fasta/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# files for recalibration
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget -P reference/broad_hg38/ -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

# files for variant annotation
time rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz reference/vep/
time tar -zxf reference/vep/homo_sapiens_vep_113_GRCh38.tar.gz -C reference/vep/
time rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/fasta/homo_sapiens/dna_index/ reference/vep/homo_sapiens/113_GRCh38/
