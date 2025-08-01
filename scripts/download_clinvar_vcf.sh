#!/bin/bash

mkdir -p reference/clinvar

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O reference/clinvar/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi -O reference/clinvar/clinvar.vcf.gz.tbi