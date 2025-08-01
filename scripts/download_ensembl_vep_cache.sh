#!/bin/bash

mkdir reference/vep

rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz reference/vep/
tar -zxf reference/vep/homo_sapiens_vep_113_GRCh38.tar.gz -C reference/vep/
rm reference/vep/homo_sapiens_vep_113_GRCh38.tar.gz
rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-113/fasta/homo_sapiens/dna_index/ reference/vep/homo_sapiens/113_GRCh38/