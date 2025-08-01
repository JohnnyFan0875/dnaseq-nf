process gatk_hard_filtering {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(raw_vcf)
    val result_dir
    path ref_fasta_dir

    output:
    tuple val(sample_id), path("${sample_id}.hardfiltered.vcf.gz"), emit: filtered_vcf

    publishDir "${result_dir}/variant_filtering", mode: 'copy'

    script:
    """
    gatk IndexFeatureFile -I ${raw_vcf}

    # Select SNPs
    gatk SelectVariants \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -V ${raw_vcf} \
        --select-type-to-include SNP \
        -O ${sample_id}.snps.vcf.gz

    # Filter SNPs with split filters
    gatk VariantFiltration \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -V ${sample_id}.snps.vcf.gz \
        --filter-name "QD_fail" --filter-expression "QD < 2.0" \
        --filter-name "FS_fail" --filter-expression "FS > 60.0" \
        --filter-name "MQ_fail" --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum_fail" --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum_fail" --filter-expression "ReadPosRankSum < -8.0" \
        -O ${sample_id}.snps.filtered.vcf.gz

    # Select Indels
    gatk SelectVariants \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -V ${raw_vcf} \
        --select-type-to-include INDEL \
        -O ${sample_id}.indels.vcf.gz

    # Filter Indels with split filters
    gatk VariantFiltration \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -V ${sample_id}.indels.vcf.gz \
        --filter-name "QD_fail" --filter-expression "QD < 2.0" \
        --filter-name "FS_fail" --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum_fail" --filter-expression "ReadPosRankSum < -20.0" \
        -O ${sample_id}.indels.filtered.vcf.gz

    # Merge filtered SNPs and Indels
    gatk MergeVcfs \
        -I ${sample_id}.snps.filtered.vcf.gz \
        -I ${sample_id}.indels.filtered.vcf.gz \
        -O ${sample_id}.hardfiltered.vcf.gz

    # Index the final VCF
    gatk IndexFeatureFile -I ${sample_id}.hardfiltered.vcf.gz
    """
}
