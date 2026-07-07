// ============================================================================
// BASEJ-DNAQC LOCAL MODULES
// ============================================================================
// Self-contained Ginkgo CNV process definitions plus the GINKO_NOPUBLISH batch
// workflow. These processes were previously included from ../../modules/ginkgo/*
// and carried a publish_dir/enable_publish/timestamp publishDir directive. They
// are inlined here WITHOUT that publishing machinery so basej-dnaqc is
// self-contained: publishing of the CNV outputs is handled by the output{} block
// in main.nf via the emitted channels, not by per-process publishDir. As a
// result the `publish_dir` param is no longer needed by this pipeline.
// ============================================================================

nextflow.enable.dsl=2

// ============================================================================
// PROCESS: BAM_TO_BED — per-sample BAM -> BED (bedtools)
// ============================================================================
process BAM_TO_BED {
    tag "${sample_name}_binsize_${bin_size}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    val(bin_size)

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

// ============================================================================
// PROCESS: GINKGO_BINUNSORT — per-sample bin unsorted reads
// ============================================================================
process GINKGO_BINUNSORT {
    tag "${sample_name}_binsize_${bin_size}"

    input:
    tuple val(sample_name), path(beds)
    path(ginkgo_binning_file)
    val(bin_size)

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

// ============================================================================
// PROCESS: GINKGO_SEGMENTATION_R — batch CNV segmentation (all samples together)
// ============================================================================
process GINKGO_SEGMENTATION_R {
    tag "binsize_${bin_size}"

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

// ============================================================================
// PROCESS: GINKGO_CNV_CALLER — batch CNV calling
// ============================================================================
process GINKGO_CNV_CALLER {
    tag "binsize_${bin_size}"

    input:
    path(SegCopy)
    val(bin_size)

    output:
    path("*.tsv"), emit: cnvs
    path("ginko_version.yml"), emit: version

    script:
    """

    CNVcaller ${SegCopy} CNV1_binsize_${bin_size}.tsv CNV2_binsize_${bin_size}.tsv

    echo Ginkgo: 0.0.2 > ginko_version.yml
    """
}

// ============================================================================
// PROCESS: GINKO_RDS_TO_FLAT — batch RDS -> flat TSVs
// ============================================================================
process GINKO_RDS_TO_FLAT {
    tag "binsize_${bin_size}"

    input:
    path (rds_file)
    val(bin_size)

    output:
    path("*.tsv"), emit: tsvs

    script:
    """

    /usr/bin/Rscript /usr/local/bin/rds_to_flat.R  ${rds_file}

    """
}

// ============================================================================
// PROCESS: GINKO_PARSE_OUTPUTS — parse Ginkgo outputs into per-sample CNV calls
// ============================================================================
process GINKO_PARSE_OUTPUTS {
    tag "binsize_${bin_size}"

    input:
    path (rds_outputs)
    path(cnvs)
    path(binref)
    val(bin_size)

    output:
    path("*.tsv"), emit: tsvs

    script:
    """

    python /scripts/divide_clouds_by_sample.py -c clouds.tsv -v1 CNV1_*.tsv -v2 CNV2_*.tsv -b ${binref}

    """
}

// ============================================================================
// PROCESS: PARSE_RDS_CNV_METRICS — batch RDS -> CNV metrics summary
// ============================================================================
process PARSE_RDS_CNV_METRICS {
    tag "PARSE_RDS_CNV_METRICS"

    input:
    path(rds_file)

    output:
    path("AllSample-GinkgoSegmentSummary.txt")

    script:
    """

    Rscript /usr/local/bin/rscript_parse_rds_cnv_metrics.R ${rds_file}


    """
}

// ============================================================================
// WORKFLOW: GINKO_NOPUBLISH
// Description: CNV calling with Ginkgo (batch processing, no publishing)
// Note: Ginkgo requires all samples together for CNV normalization
// ============================================================================
workflow GINKO_NOPUBLISH {
    take:
        ch_bam
        ch_bin_size
        ch_binref
        ch_gcref
        ch_boundsref_file
        ch_min_ploidy
        ch_max_ploidy
        ch_min_bin_width
        ch_is_haplotype

    main:
        // Per-sample: Convert BAM to BED
        BAM_TO_BED(
            ch_bam,
            ch_bin_size
        )

        // Per-sample: Bin unsorted reads
        GINKGO_BINUNSORT(
            BAM_TO_BED.out.bed_only,
            ch_binref.collect(),
            ch_bin_size
        )

        // Collect all mapped files for batch processing
        ch_mapped_files = GINKGO_BINUNSORT.out.map { it -> it.last() }.collect()

        // Save mapped files channel for raw_counts output
        ch_raw_counts = ch_mapped_files

        // Batch: CNV segmentation (requires all samples together)
        GINKGO_SEGMENTATION_R(
            ch_mapped_files,
            ch_binref.collect(),
            ch_gcref.collect(),
            ch_boundsref_file.collect(),
            ch_min_ploidy,
            ch_max_ploidy,
            ch_min_bin_width,
            ch_bin_size,
            ch_is_haplotype
        )

        // Batch: Call CNVs
        GINKGO_CNV_CALLER(
            GINKGO_SEGMENTATION_R.out.segcopy,
            ch_bin_size
        )

        // Batch: Convert RDS to flat files
        GINKO_RDS_TO_FLAT(
            GINKGO_SEGMENTATION_R.out.RDS,
            ch_bin_size
        )

        // Batch: Parse Ginkgo outputs
        GINKO_PARSE_OUTPUTS(
            GINKO_RDS_TO_FLAT.out.tsvs,
            GINKGO_CNV_CALLER.out.cnvs,
            ch_binref.collect(),
            ch_bin_size
        )

        // Batch: Parse RDS for CNV metrics summary
        PARSE_RDS_CNV_METRICS(
            GINKGO_SEGMENTATION_R.out.RDS
        )

    emit:
        // AllSample-GinkgoSegmentSummary.txt (one row per sample)
        metrics = PARSE_RDS_CNV_METRICS.out
        // Per-sample CNV calls
        cnvs = GINKO_PARSE_OUTPUTS.out.tsvs
        // CNV graph
        graph = GINKGO_SEGMENTATION_R.out.jpeg
        // RDS file (batch)
        rds = GINKGO_SEGMENTATION_R.out.RDS
        // SegCopy file (batch - samples as columns)
        segcopy = GINKGO_SEGMENTATION_R.out.segcopy
        // Raw counts (mapped files for parquet generation)
        raw_counts = ch_raw_counts
        // Merged raw counts file (binUnsorted.outdata) - matches SegCopy structure
        raw_counts_merged = GINKGO_SEGMENTATION_R.out.raw_counts_merged
        // Version info
        ginkgo_version = GINKGO_CNV_CALLER.out.version
        bedtools_version = BAM_TO_BED.out.version
}
