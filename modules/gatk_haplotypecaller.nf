process gatk_haplotypecaller {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_bam), path(recal_bai)
    val result_dir
    path ref_fasta_dir

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf_variant_calls
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf_variant_calls

    publishDir "${result_dir}/variant_calls", mode: 'copy'

    script:
    """
    # Generate gVCF
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -I ${recal_bam} \
        -O ${sample_id}.g.vcf.gz \
        -ERC GVCF

    # Generate VCF
    gatk --java-options "-Xmx8g" GenotypeGVCFs \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -V ${sample_id}.g.vcf.gz \
        -O ${sample_id}.vcf.gz    
    """
}
