nextflow.enable.dsl=2
params.timestamp = ""

process GINKO_PARSE_OUTPUTS {
    tag "binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    path (rds_outputs)
    path(cnvs)
    path(binref)
    val(bin_size)
    val(publish_dir)
    val(enable_publish)
    

    output:
    path("*.tsv"), emit: tsvs
    
    script:
    """
        
    python /scripts/divide_clouds_by_sample.py -c clouds.tsv -v1 CNV1_*.tsv -v2 CNV2_*.tsv -b ${binref}
        

    """
}

