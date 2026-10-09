nextflow.enable.dsl=2
params.timestamp = ""

process GINKGO_BINUNSORT {
    tag "${sample_name}_binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    tuple val(sample_name), path(beds)
    path(ginkgo_binning_file)
    val(bin_size)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("*.mapped"), emit: binunsorted_data

    
    script:
    """
    #!/bin/bash
    NB_BINS=\$(wc -l < ${ginkgo_binning_file})
    name=\$(echo ${beds}  | sed 's|.bed||')
    binUnsorted ${ginkgo_binning_file} \${NB_BINS} ${beds} \${name} \${name}.mapped

    """
}

workflow GINKGO_BINUNSORT_WF {
    take:
        ch_beds
        ch_ginkgo_binning_file
        ch_bin_size
        ch_publish_dir
        ch_enable_publish
        
    main:
        GINKGO_BINUNSORT(
                            ch_beds,
                            ch_ginkgo_binning_file,
                            ch_bin_size,
                            ch_publish_dir,
                            ch_enable_publish
                        )
    emit:
        binunsorted_data = GINKGO_BINUNSORT.out.binunsorted_data
}

workflow{
    
    ch_beds = Channel.fromFilePairs(params.beds + "/*.bed")
    ch_ginkgo_binning_file = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_76_bwa")
    
    GINKGO_BINUNSORT_WF (
                            ch_beds,
                            ch_ginkgo_binning_file.collect(),
                            params.bin_size,
                            params.publish_dir,
                            params.enable_publish
                        )
}
