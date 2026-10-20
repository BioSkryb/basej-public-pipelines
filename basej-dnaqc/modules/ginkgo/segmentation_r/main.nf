nextflow.enable.dsl=2
params.timestamp = ""

process GINKGO_SEGMENTATION_R {
    tag "binsize_${bin_size}"
    publishDir  "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginkgo", enabled:"$enable_publish"
    
    input:
    path(binunsorted_file)
    path(binref_file)
    path(gcref_file)
    path(boundsref_file)
    val(min_ploidy)
    val(max_ploidy)
    val(min_bin_width)
    val(bin_size)
    val(is_haplotype)
    val(publish_dir)
    val(enable_publish)


    output:
    path("*ginkgo_res.binsize_${bin_size}.RDS"), emit: RDS
    path("cnv_plots_binsize_${bin_size}.tar.gz"), emit: jpeg
    path("*SegCopy.binsize_${bin_size}.tsv"), emit: segcopy
    path("binUnsorted.${bin_size}.outdata"), emit: raw_counts_merged

    script:
    """
    paste *.mapped > binUnsorted.${bin_size}.outdata
    if [ ${is_haplotype} == 0 ]; then
    
        echo "Other";    
        
        /usr/bin/Rscript /usr/local/bin/cnv_ginkgo.R binUnsorted.${bin_size}.outdata ${binref_file} ${gcref_file} ${boundsref_file} ${min_ploidy} ${max_ploidy} ${min_bin_width} 0
        
        
    else
    
        echo "Provided ploidy";
    
    cat binUnsorted.${bin_size}.outdata | head -n1 | tr '\\t' '\\n' | while read line; 
    do      
        
        echo -e "\${line}\t${is_haplotype}" 
        
    done > ploidy.txt
    
    /usr/bin/Rscript /usr/local/bin/cnv_ginkgo.R binUnsorted.${bin_size}.outdata ${binref_file} ${gcref_file} ${boundsref_file} ${min_ploidy} ${max_ploidy} ${min_bin_width} 1
    
    fi
    
    mv ginkgo_res.RDS ginkgo_res.binsize_${bin_size}.RDS
    mv SegCopy SegCopy.binsize_${bin_size}.tsv
    
    ls *.jpeg | while read file;
    do 
    
        mv "\${file}" cnv_binsize_${bin_size}_"\${file}"
        
    done
    
    # Compress all JPEG files into a tar.gz archive
    tar -czf cnv_plots_binsize_${bin_size}.tar.gz cnv_binsize_${bin_size}_*.jpeg
    
    """
}

workflow GINKGO_SEGMENTATION_R_WF {
    take:
        ch_binunsorted
        ch_binref
        ch_gcref
        ch_boundsref_file
        ch_min_ploidy
        ch_max_ploidy
        ch_min_bin_width
        ch_bin_size
        ch_is_haplotype
        ch_publish_dir
        ch_enable_publish
        
    main:
        GINKGO_SEGMENTATION_R(
                               ch_binunsorted,
                               ch_binref,
                               ch_gcref,
                               ch_boundsref_file,
                               ch_min_ploidy,
                               ch_max_ploidy,
                               ch_min_bin_width,
                               ch_bin_size,
                               ch_is_haplotype,
                               ch_publish_dir,
                               ch_enable_publish 
                             )
        
    emit:
        RDS = GINKGO_SEGMENTATION_R.out.RDS
        jpeg = GINKGO_SEGMENTATION_R.out.jpeg
        segcopy = GINKGO_SEGMENTATION_R.out.segcopy
        raw_counts_merged = GINKGO_SEGMENTATION_R.out.raw_counts_merged
}

workflow{
    
    ch_binunsorted = Channel.fromFilePairs(params.binunsorted + "/*.outdata")
    ch_binref = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_76_bwa")
    ch_gcref = Channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_76_bwa")
    ch_boundsref_file = Channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_76_bwa")
    
    GINKGO_SEGMENTATION_R_WF (
                               ch_binunsorted,
                               ch_binref.collect(),
                               ch_gcref.collect(),
                               ch_boundsref_file.collect(),
                               params.min_ploidy,
                               params.max_ploidy,
                               params.min_bin_width,
                               params.bin_size,
                               params.is_haplotype,
                               params.publish_dir,
                               params.enable_publish 
                             )
}
