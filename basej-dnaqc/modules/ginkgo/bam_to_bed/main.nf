nextflow.enable.dsl=2
params.timestamp = ""

process BAM_TO_BED {
    tag "${sample_name}_binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    val(bin_size)
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("${sample_name}.bed"), emit: bed_only
    path("bedtools_version.yml"), emit: version

    
    script:
    """
    bedtools bamtobed -i ${bam[0]} > ${sample_name}.bed
    
    export BEDTOOLS_VER=\$(bedtools --version 2>&1 | sed -e "s/bedtools //g")
    echo bedtools: \$BEDTOOLS_VER > bedtools_version.yml
    """
}

workflow BAM_TO_BED_WF {
    take:
       ch_bam
       ch_bin_size
       ch_publish_dir
       ch_enable_publish
        
    main:
        BAM_TO_BED(
                    ch_bam,
                    ch_bin_size,
                    ch_publish_dir,
                    ch_enable_publish
                  )
    emit:
        bed_only = BAM_TO_BED.out.bed_only
        version = BAM_TO_BED.out.version
}

workflow {
    ch_bam_raw = Channel.fromFilePairs(params.bam_dir + "/*{.bam,bam.bai}", size: -1)
    ch_bam_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_bam }
              
    BAM_TO_BED_WF(
                    ch_bam,
                    params.bin_size,
                    params.publish_dir,
                    params.enable_publish
                 )
}
