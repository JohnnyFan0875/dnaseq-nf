process mark_duplicates {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    val result_dir

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bai"), path("${sample_id}.metrics.txt"), emit: dedup_bams

    publishDir { "${result_dir}/dedup" }, mode: 'rellink'

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \
        I=${bam} \
        O=${sample_id}.dedup.bam \
        M=${sample_id}.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT
    """
}
