nextflow.enable.dsl=2
params.timestamp = ""

process GINKO_RDS_TO_FLAT {
    tag "binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    path (rds_file)
    val(bin_size)
    val(publish_dir)
    val(enable_publish)
    

    output:
    path("*.tsv"), emit: tsvs
    
    script:
    """
        
    /usr/bin/Rscript /usr/local/bin/rds_to_flat.R  ${rds_file}
        

    """
}

