process vep_annotation {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(filtered_vcf)
    val result_dir
    path ref_fasta_dir
    path vep_cache_dir
    path clinvar_dir

    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"), emit: vep_annotated_vcf

    publishDir "${result_dir}/variant_annotation", mode: 'copy'

    script:
    """
    vep \
      --input_file ${filtered_vcf} \
      --output_file ${sample_id}.vep.vcf.gz \
      --vcf \
      --compress_output bgzip \
      --dir_cache ${vep_cache_dir} \
      --fasta ${ref_fasta_dir}/Homo_sapiens_assembly38.fasta \
      --assembly GRCh38 \
      --custom ${clinvar_dir}/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
      --force_overwrite \
      --offline

    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
}
