process gatk_haplotypecaller {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_bam), path(recal_bai)
    val result_dir
    path ref_fasta_dir

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), optional: true, emit: gvcf_variant_calls
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf_variant_calls

    publishDir "${result_dir}/variant_calls", mode: 'copy'

    script:
    def intervals = params.test_intervals ? "-L ${params.test_intervals}" : ""
    def erc_flag  = params.test_skip_gvcf ? "" : "-ERC GVCF"

    """
    # Generate variant file
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -I ${recal_bam} \
        ${intervals} \
        ${erc_flag} \
        -O ${sample_id}.${params.test_skip_gvcf ? 'vcf' : 'g.vcf'}.gz

    # Generate VCF if gVCF is not skipped
    ${params.test_skip_gvcf ? "" : """
        gatk --java-options "-Xmx8g" GenotypeGVCFs \
            -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
            -V ${sample_id}.g.vcf.gz \
            -O ${sample_id}.vcf.gz
        """
    }
    """
}
