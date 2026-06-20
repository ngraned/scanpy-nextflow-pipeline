process scanpy_analysis {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
        path matrix_dir
    
    output:
        path "adata_processed.h5ad"
        path "umap.png"
        path "report.txt"
    
    script:
    """
    python ${projectDir}/scanpy_analysis.py ${matrix_dir}
    """
}

workflow {
    if (!params.input_mtx) {
        error("--input_mtx is required")
    }
    scanpy_analysis(params.input_mtx)
}
