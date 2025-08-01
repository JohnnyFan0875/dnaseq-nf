nextflow.enable.dsl=2

// Include modules
include { fastqc as fastqc_pre } from './modules/fastqc.nf'
include { fastqc as fastqc_post } from './modules/fastqc.nf'
include { fastp }  from './modules/fastp.nf'
include { multiqc } from './modules/multiqc.nf'
include { align_sort_reads } from './modules/align_sort_reads.nf'
include { mark_duplicates } from './modules/markduplicates.nf'
include { base_recalibration } from './modules/base_recalibration.nf'
include { gatk_haplotypecaller } from './modules/gatk_haplotypecaller.nf'
include { gatk_hard_filtering } from './modules/gatk_hard_filtering.nf'
include { vep_annotation } from './modules/vep_annotation.nf'

// Define parameter
if (!params.project_name) {
    error "Please provide a project name using --project_name"
}

params.project_dir = "./data/${params.project_name}"
params.samplesheet = "${params.project_dir}/metadata/sample_metadata.csv"
params.result_dir = "${params.project_dir}/results/"
params.ref_dir = "./reference/"
params.vep_cache_dir = "${params.ref_dir}/vep_cache/"

workflow {
    samples = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            tuple(row.sample_id, 
                  file("${params.project_dir}/raw/${row.fastq1}"), 
                  file("${params.project_dir}/raw/${row.fastq2}"))
        }

    // Run FastQC (pre-trim)
    fastqc_pre_out = fastqc_pre(samples, Channel.value("pre"), Channel.value(params.result_dir))

    // Run Fastp
    trimmed_out = fastp(samples, Channel.value(params.result_dir))

    // Run FastQC (post-trim)
    fastqc_post_out = fastqc_post(trimmed_out.trimmed_reads, Channel.value("post"), Channel.value(params.result_dir))

    // multiqc (without post-trimming fastqc)
    all_qc_reports = fastqc_pre_out.fastqc_reports
        .mix(trimmed_out.html_reports)
        .mix(trimmed_out.json_reports)
        .collect()
 
    multiqc(all_qc_reports, Channel.value(params.result_dir))
 
    // Align and sort reads
    aligned_bams = align_sort_reads(trimmed_out.trimmed_reads, 
                                    Channel.value(params.result_dir),
                                    Channel.fromPath("${params.ref_dir}/fasta/"))

    // Mark duplicates
    dedup_bams = mark_duplicates(aligned_bams.aligned_reads,
                                 Channel.value(params.result_dir))

    // Base recalibration
    recal_bams = base_recalibration(dedup_bams.dedup_bams,
                                    Channel.value(params.result_dir),
                                    Channel.fromPath("${params.ref_dir}/fasta/"),
                                    Channel.fromPath("${params.ref_dir}/broad_hg38/"))

    // Variant calling
    vcf_files = gatk_haplotypecaller(recal_bams.recal_bams,
                                     Channel.value(params.result_dir),
                                     Channel.fromPath("${params.ref_dir}/fasta/"))

    // Hard filtering
    filtered_vcfs = gatk_hard_filtering(vcf_files.vcf_variant_calls,
                                        Channel.value(params.result_dir),
                                        Channel.fromPath("${params.ref_dir}/fasta/"))

    // VEP annotation
    vep_annotated_vcfs = vep_annotation(filtered_vcfs.filtered_vcf,
                                        Channel.value(params.result_dir),
                                        Channel.fromPath("${params.ref_dir}/fasta/"),
                                        Channel.fromPath("${params.ref_dir}/vep/"),
                                        Channel.fromPath("${params.ref_dir}/clinvar/"))
}