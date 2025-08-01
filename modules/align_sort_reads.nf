process align_sort_reads {

    tag "$sample_id"

    cpus 8

    input:
    tuple val(sample_id), path(read1), path(read2)
    val result_dir
    path ref_fasta_dir

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: aligned_reads

    publishDir { "${result_dir}/bwa_align" }, mode: 'rellink'

    script:
    """
    id=\$(zcat "${read1}" | head -n 1 | cut -f 3-4 -d':' | sed 's/@//') # flowcell ID + lane
    bwa mem -M -R "@RG\\tID:\${id}\\tSM:${sample_id}\\tPL:ILLUMINA" -t ${task.cpus} "${ref_fasta_dir}/Homo_sapiens_assembly38.fasta" "${read1}" "${read2}" | \
    samtools sort -@${task.cpus} -o "${sample_id}.sorted.bam" -
    samtools index "${sample_id}.sorted.bam"
    """
}
