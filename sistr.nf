#!/usr/bin/env nextflow

// Input parameters
params.input_s3 = "s3://your-bucket-name/validationset"
params.outdir = "s3://your-bucket-name/sistr_results"

workflow {
    def contig_files = channel.fromPath("${params.input_s3}/*.contigs.fa")
        .map { file -> [file.baseName.replaceAll('.contigs.fa', ''), file] }

    runSISTR(contig_files)
    combineResults(runSISTR.out.collect())
}

process runSISTR {
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample_id), path(contig_file)
    
    output:
    path("${sample_id}_sistr_results.txt")
    
    script:
    """
    sistr --qc-filter off -i ${contig_file} -o ${sample_id}_sistr_results.txt
    """
}

process combineResults {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path(all_results)
    
    output:
    path("combined_sistr_results.tab")
    
    script:
    """
    cat ${all_results} > combined_sistr_results.tab
    """
}
