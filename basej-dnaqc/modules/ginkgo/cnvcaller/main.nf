nextflow.enable.dsl=2
params.timestamp = ""

process GINKGO_CNV_CALLER {
    tag "binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    path(SegCopy)
    val(bin_size)
    val(publish_dir)
    val(enable_publish)

    output:
    path("*.tsv"), emit: cnvs
    path("ginko_version.yml"), emit: version

    
    script:
    """
    
    CNVcaller ${SegCopy} CNV1_binsize_${bin_size}.tsv CNV2_binsize_${bin_size}.tsv
    
    echo Ginkgo: 0.0.2 > ginko_version.yml
    """
}

workflow GINKGO_CNV_CALLER_WF {
    take:
        ch_SegCopy
        ch_bin_size
        ch_publish_dir
        ch_enable_publish
        
    main:
        GINKGO_CNV_CALLER(
                            ch_SegCopy,
                            ch_bin_size,
                            ch_publish_dir,
                            ch_enable_publish
                         )
        
    emit:
        cnvs = GINKGO_CNV_CALLER.out.cnvs
        version = GINKGO_CNV_CALLER.out.version
}

workflow {
    
    ch_SegCopy  = Channel.fromFilePairs(params.SegCopy + "/*.tsv")
    GINKGO_CNV_CALLER_WF (
                            ch_SegCopy,
                            params.bin_size,
                            params.publish_dir,
                            params.enable_publish
                          )
}
