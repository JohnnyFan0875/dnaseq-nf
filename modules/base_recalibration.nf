process base_recalibration {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai), path(metrics)
    val result_dir
    path ref_fasta_dir
    path ref_broad_dir

    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bai"), emit: recal_bams

    publishDir "${result_dir}/recalibrated", mode: 'rellink'

    script:
    """
    gatk BaseRecalibrator \
        -I ${dedup_bam} \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        --known-sites ${ref_broad_dir}/Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ${ref_broad_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -O ${sample_id}.recal.table

    gatk ApplyBQSR \
        -R ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
        -I ${dedup_bam} \
        --bqsr-recal-file ${sample_id}.recal.table \
        -O ${sample_id}.recal.bam

    gatk BuildBamIndex -I ${sample_id}.recal.bam
    """
}
