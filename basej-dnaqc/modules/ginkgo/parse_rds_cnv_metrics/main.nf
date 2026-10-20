nextflow.enable.dsl=2
params.timestamp = ""

process PARSE_RDS_CNV_METRICS {
    tag "PARSE_RDS_CNV_METRICS"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    path(rds_file)
    val(publish_dir)
    val(enable_publish)

    output:
    path("AllSample-GinkgoSegmentSummary.txt")

    script:
    """
    
    Rscript /usr/local/bin/rscript_parse_rds_cnv_metrics.R ${rds_file}
   
    
    """
}
