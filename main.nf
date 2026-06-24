process scanpy_analysis {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
        path scanpy_script
        path matrix_dir
    
    output:
        path "adata_processed.h5ad"
        path "umap.png"
        path "report.txt"
    
    script:
    """
    python ${scanpy_script} ${matrix_dir}
    """
}

workflow {
    if (!params.input_mtx) {
        error("--input_mtx is required")
    }
    script_ch = Channel.fromPath("${projectDir}/scanpy_analysis.py")
    scanpy_analysis(script_ch, params.input_mtx)
}
