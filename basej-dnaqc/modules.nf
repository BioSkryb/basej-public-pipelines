// ============================================================================
// BASEJ-METHYLQC LOCAL MODULES
// ============================================================================
// Cloned from basej-dnaqc - Ginkgo CNV workflow for conditional CNV calling
// ============================================================================

nextflow.enable.dsl=2

// ============================================================================
// INCLUDE: GINKGO MODULE DEPENDENCIES
// ============================================================================
include { BAM_TO_BED } from './modules/ginkgo/bam_to_bed/main.nf'
include { GINKGO_BINUNSORT } from './modules/ginkgo/binunsort/main.nf'
include { GINKGO_SEGMENTATION_R } from './modules/ginkgo/segmentation_r/main.nf'
include { GINKGO_CNV_CALLER } from './modules/ginkgo/cnvcaller/main.nf'
include { GINKO_RDS_TO_FLAT } from './modules/ginkgo/rds_to_flat/main.nf'
include { PARSE_RDS_CNV_METRICS } from './modules/ginkgo/parse_rds_cnv_metrics/main.nf'
include { GINKO_PARSE_OUTPUTS } from './modules/ginkgo/parse_ginko_outputs/main.nf'

// ============================================================================
// WORKFLOW: GINKO_NOPUBLISH
// Description: CNV calling with Ginkgo (batch processing, no publishing)
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
        BAM_TO_BED(ch_bam, ch_bin_size, ".", false)
        GINKGO_BINUNSORT(BAM_TO_BED.out.bed_only, ch_binref.collect(), ch_bin_size, ".", false)
        ch_mapped_files = GINKGO_BINUNSORT.out.map { it -> it.last() }.collect()
        ch_raw_counts = ch_mapped_files

        GINKGO_SEGMENTATION_R(
            ch_mapped_files,
            ch_binref.collect(),
            ch_gcref.collect(),
            ch_boundsref_file.collect(),
            ch_min_ploidy,
            ch_max_ploidy,
            ch_min_bin_width,
            ch_bin_size,
            ch_is_haplotype,
            ".",
            false
        )
        GINKGO_CNV_CALLER(GINKGO_SEGMENTATION_R.out.segcopy, ch_bin_size, ".", false)
        GINKO_RDS_TO_FLAT(GINKGO_SEGMENTATION_R.out.RDS, ch_bin_size, ".", false)
        GINKO_PARSE_OUTPUTS(
            GINKO_RDS_TO_FLAT.out.tsvs,
            GINKGO_CNV_CALLER.out.cnvs,
            ch_binref.collect(),
            ch_bin_size,
            ".",
            false
        )
        PARSE_RDS_CNV_METRICS(GINKGO_SEGMENTATION_R.out.RDS, ".", false)

    emit:
        metrics = PARSE_RDS_CNV_METRICS.out
        cnvs = GINKO_PARSE_OUTPUTS.out.tsvs
        graph = GINKGO_SEGMENTATION_R.out.jpeg
        rds = GINKGO_SEGMENTATION_R.out.RDS
        segcopy = GINKGO_SEGMENTATION_R.out.segcopy
        raw_counts = ch_raw_counts
        raw_counts_merged = GINKGO_SEGMENTATION_R.out.raw_counts_merged
        ginkgo_version = GINKGO_CNV_CALLER.out.version
        bedtools_version = BAM_TO_BED.out.version
}
