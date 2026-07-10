nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// MUTSIGNATURES_FROM_MATRICES  (named subworkflow, basej style)
//
// COSMIC mutational signature analysis directly from a pre-binarized variant ×
// sample matrix — the "binary_matrix" input mode (Option D). No VCFs are read.
//
// Input (REQUIRED): binary_matrix — pre-binarized variant × sample TSV with 0/1
//   values; row keys {chrom}_{pos}_{ref}_{alt} (pos = rightmost numeric token);
//   columns = sample names.
//
// Pipeline:
//   EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX (per-sample variant tables)
//     → SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES (SBS96, DBS78, ID83)
//     → SIGPROFILER_ASSIGNMENT_{SBS,DBS,ID} (cohort cosmic_fit per mutation type)
//     → MERGE_ACTIVITIES_{SBS,DBS,ID}
//     → COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE (per sample; combined SBS+DBS+ID dirs)
//     → PLOT_BARGRAPHS_{SBS,DBS,ID} (stacked-bar PNG/PDF + HTML per type)
//     → GATHER_MUTSIG_ARTIFACTS
//
// Optional (toggled by take flags): MUSICAL / SIGDYN / MUSICATK, and the comparative
// plot (only when both MuSiCaL and SigDyn are enabled).
//
// DBS and ID branches are conditional: they only execute when the corresponding
// matrix from SigProfilerMatrixGenerator is non-empty. SBS always executes.
// ─────────────────────────────────────────────────────────────────────────────

include { EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX          } from '../modules.nf'
include { SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES } from '../modules.nf'
include { COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE                      } from '../modules.nf'
include { GATHER_MUTSIG_ARTIFACTS                                         } from '../modules.nf'

// Triple-aliased includes for SBS, DBS, and ID mutation types
include { SIGPROFILER_ASSIGNMENT_FROM_MATRIX as SIGPROFILER_ASSIGNMENT_SBS } from '../modules.nf'
include { SIGPROFILER_ASSIGNMENT_FROM_MATRIX as SIGPROFILER_ASSIGNMENT_DBS } from '../modules.nf'
include { SIGPROFILER_ASSIGNMENT_FROM_MATRIX as SIGPROFILER_ASSIGNMENT_ID  } from '../modules.nf'

include { MERGE_SIGNATURE_ACTIVITIES as MERGE_ACTIVITIES_SBS } from '../modules.nf'
include { MERGE_SIGNATURE_ACTIVITIES as MERGE_ACTIVITIES_DBS } from '../modules.nf'
include { MERGE_SIGNATURE_ACTIVITIES as MERGE_ACTIVITIES_ID  } from '../modules.nf'

include { PLOT_MATRIX_SIGNATURE_BARGRAPHS as PLOT_BARGRAPHS_SBS } from '../modules.nf'
include { PLOT_MATRIX_SIGNATURE_BARGRAPHS as PLOT_BARGRAPHS_DBS } from '../modules.nf'
include { PLOT_MATRIX_SIGNATURE_BARGRAPHS as PLOT_BARGRAPHS_ID  } from '../modules.nf'

// Additional signature analysis pipelines (toggled by musical_enabled / sigdyn_enabled / musicatk_enabled)
include { MUSICAL_SIGNATURE_ANALYSIS            } from '../modules.nf'
include { SIGDYN_SIGNATURE_ANALYSIS             } from '../modules.nf'
include { MUSICATK_SIGNATURE_ANALYSIS           } from '../modules.nf'
include { PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES } from '../modules.nf'

// Lightweight process to check if a matrix file exists and is non-empty.
// Runs inside the execution environment where matrix_generator_out is staged,
// avoiding issues with S3 path resolution in channel .map operators.
// Emits the matrix directory path as a string via stdout only when the target
// matrix file has content (>1 line, i.e., header + data).
// Container + resources are configured in nextflow.config under withName.
process CHECK_MATRIX_NON_EMPTY_DBS {
    tag "DBS"

    input:
    path(matrix_dir)
    val(cohort_id)
    val(scope_suffix)

    output:
    stdout emit: result

    script:
    """
    f="${matrix_dir}/DBS/${cohort_id}.DBS78.${scope_suffix}"
    if [ -f "\$f" ] && [ -s "\$f" ]; then
        lines=\$(wc -l < "\$f")
        if [ "\$lines" -gt 1 ]; then
            printf "HAS_DATA"
        else
            printf "EMPTY"
        fi
    else
        printf "EMPTY"
    fi
    """
}

process CHECK_MATRIX_NON_EMPTY_ID {
    tag "ID"

    input:
    path(matrix_dir)
    val(cohort_id)
    val(scope_suffix)

    output:
    stdout emit: result

    script:
    """
    f="${matrix_dir}/ID/${cohort_id}.ID83.${scope_suffix}"
    if [ -f "\$f" ] && [ -s "\$f" ]; then
        lines=\$(wc -l < "\$f")
        if [ "\$lines" -gt 1 ]; then
            printf "HAS_DATA"
        else
            printf "EMPTY"
        fi
    else
        printf "EMPTY"
    fi
    """
}

workflow MUTSIGNATURES_FROM_MATRICES {

    take:
    ch_binary                              // channel: pre-binarized variant × sample matrix
    cohort_id                              // val
    sig_genome_build                       // val
    sig_context_type                       // val
    sig_exome                              // val
    sig_export_probabilities               // val
    sig_export_probabilities_per_mutation  // val
    musical_enabled                        // val
    sigdyn_enabled                         // val
    musicatk_enabled                       // val
    musicatk_k_denovo                      // val

    main:

    // ── Step 1: per-sample variant tables from the binary matrix ──────────────
    EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX (
        ch_binary,
        cohort_id,
        sig_genome_build
    )

    def ch_per_sample_tsv = EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX.out.per_sample_tsvs
        .flatten()
        .map { tsv ->
            def sample_name = tsv.name.replaceAll('_mutsig_somatic_variant_table\\.tsv$', '')
            tuple( sample_name, tsv )
        }

    // ── Step 2: SigProfilerMatrixGenerator (SBS96, DBS78, ID83) ───────────────
    SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES (
        ch_per_sample_tsv.map { _s, tsv -> tsv }.collect(),
        cohort_id,
        sig_genome_build,
        sig_exome
    )

    // The matrix generator produces SBS, DBS, and ID subdirectories. SBS always
    // has data; DBS and ID may be empty. Convert to a value channel (.first()) so
    // it can be consumed by multiple downstream processes/operators without draining.
    ch_matrix_out = SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES.out.matrix_output.first()

    def matrixScopeSuffix = sig_exome.toString() == 'true' ? 'exome' : 'all'

    // ── Empty-matrix gating for DBS / ID branches (SBS always runs) ───────────
    CHECK_MATRIX_NON_EMPTY_DBS (
        ch_matrix_out,
        cohort_id,
        matrixScopeSuffix
    )
    ch_dbs_matrix = ch_matrix_out
        .combine(CHECK_MATRIX_NON_EMPTY_DBS.out.result)
        .filter { _dir, result -> result.trim() == 'HAS_DATA' }
        .map { dir, _result -> dir }

    CHECK_MATRIX_NON_EMPTY_ID (
        ch_matrix_out,
        cohort_id,
        matrixScopeSuffix
    )
    ch_id_matrix = ch_matrix_out
        .combine(CHECK_MATRIX_NON_EMPTY_ID.out.result)
        .filter { _dir, result -> result.trim() == 'HAS_DATA' }
        .map { dir, _result -> dir }

    // ── Step 3: SBS signature assignment (always runs) ───────────────────────
    SIGPROFILER_ASSIGNMENT_SBS (
        ch_matrix_out,
        "SBS",
        cohort_id,
        sig_genome_build,
        sig_context_type,
        sig_exome,
        sig_export_probabilities,
        sig_export_probabilities_per_mutation
    )

    MERGE_ACTIVITIES_SBS (
        SIGPROFILER_ASSIGNMENT_SBS.out.assignment_cohort_dir.collect(),
        "SBS"
    )

    // ── Step 4: DBS signature assignment (conditional on non-empty DBS matrix) ─
    SIGPROFILER_ASSIGNMENT_DBS (
        ch_dbs_matrix,
        "DBS",
        cohort_id,
        sig_genome_build,
        sig_context_type,
        sig_exome,
        sig_export_probabilities,
        sig_export_probabilities_per_mutation
    )

    MERGE_ACTIVITIES_DBS (
        SIGPROFILER_ASSIGNMENT_DBS.out.assignment_cohort_dir.collect(),
        "DBS"
    )

    // ── Step 5: ID signature assignment (conditional on non-empty ID matrix) ──
    SIGPROFILER_ASSIGNMENT_ID (
        ch_id_matrix,
        "ID",
        cohort_id,
        sig_genome_build,
        sig_context_type,
        sig_exome,
        sig_export_probabilities,
        sig_export_probabilities_per_mutation
    )

    MERGE_ACTIVITIES_ID (
        SIGPROFILER_ASSIGNMENT_ID.out.assignment_cohort_dir.collect(),
        "ID"
    )

    // ── Step 6: Per-sample signature coverage ────────────────────────────────
    // Combine all available assignment directories. SBS is always present; DBS/ID
    // only run when their matrices are non-empty.
    def ch_all_assignment_dirs = SIGPROFILER_ASSIGNMENT_SBS.out.assignment_cohort_dir
        .mix(
            SIGPROFILER_ASSIGNMENT_DBS.out.assignment_cohort_dir,
            SIGPROFILER_ASSIGNMENT_ID.out.assignment_cohort_dir
        )
        .collect()

    // Wrap the collected dir list so combine keeps it as a single tuple element
    // (path(assignment_dirs)) instead of spreading it into separate positional args.
    def ch_mat_coverage_input = ch_per_sample_tsv
        .combine(ch_all_assignment_dirs.map { dirs -> [ dirs ] })

    COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE ( ch_mat_coverage_input )

    def ch_cov_tsv_list = COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE.out.coverage_tsv
        .map { _sample_name, tsv -> tsv }.collect()

    // ── Step 7: Signature bar graphs per mutation type ────────────────────────
    // Wrap the collected coverage TSV list so combine keeps it as a single tuple
    // element (path(mutsig_coverage_tsvs)) rather than spreading it.
    def ch_cov_tsv_wrapped = ch_cov_tsv_list.map { tsvs -> [ tsvs ] }

    ch_sbs_bargraph_input = MERGE_ACTIVITIES_SBS.out.merged_activities.combine(ch_cov_tsv_wrapped)
    PLOT_BARGRAPHS_SBS ( ch_sbs_bargraph_input, "SBS" )

    ch_dbs_bargraph_input = MERGE_ACTIVITIES_DBS.out.merged_activities.combine(ch_cov_tsv_wrapped)
    PLOT_BARGRAPHS_DBS ( ch_dbs_bargraph_input, "DBS" )

    ch_id_bargraph_input = MERGE_ACTIVITIES_ID.out.merged_activities.combine(ch_cov_tsv_wrapped)
    PLOT_BARGRAPHS_ID ( ch_id_bargraph_input, "ID" )

    // ── Step 8: Additional signature analysis pipelines ───────────────────────
    def ch_musical_activities  = null
    def ch_sigdyn_activities   = null
    def ch_musicatk_pub        = channel.empty()

    // MuSiCaL — COSMIC v3.2 SBS refitting (likelihood-bidirectional NNLS)
    if ( musical_enabled.toString() == 'true' ) {
        MUSICAL_SIGNATURE_ANALYSIS (
            ch_matrix_out,
            cohort_id,
            sig_genome_build,
            sig_exome
        )
        ch_musical_activities = MUSICAL_SIGNATURE_ANALYSIS.out.musical_activities
    }

    // SigDyn — COSMIC v3.2 strict best-subset fitting + dynamics statistics
    if ( sigdyn_enabled.toString() == 'true' ) {
        SIGDYN_SIGNATURE_ANALYSIS (
            ch_matrix_out,
            cohort_id,
            sig_genome_build,
            sig_exome
        )
        ch_sigdyn_activities = SIGDYN_SIGNATURE_ANALYSIS.out.sigdyn_activities
    }

    // musicatk — de novo NMF/LDA signature discovery from per-sample variant TSVs
    if ( musicatk_enabled.toString() == 'true' ) {
        def ch_musicatk_tsvs = ch_per_sample_tsv
            .map { _s, tsv -> tsv }
            .collect()

        MUSICATK_SIGNATURE_ANALYSIS (
            ch_musicatk_tsvs,
            cohort_id,
            sig_genome_build,
            sig_exome,
            musicatk_k_denovo
        )
        ch_musicatk_pub = MUSICATK_SIGNATURE_ANALYSIS.out.musicatk_activities
            .mix( MUSICATK_SIGNATURE_ANALYSIS.out.musicatk_signatures )
    }

    // ── Step 9: Comparative faceted signature activity plot ──────────────────
    // Runs when both MuSiCaL and SigDyn are enabled (3-panel: SPA, MuSiCaL, SigDyn).
    def ch_comparative_figs = channel.empty()
    if ( ch_musical_activities && ch_sigdyn_activities ) {
        PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES (
            MERGE_ACTIVITIES_SBS.out.merged_activities,
            ch_musical_activities,
            ch_sigdyn_activities,
            cohort_id
        )
        ch_comparative_figs = PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES.out.comparative_png
            .mix(
                PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES.out.comparative_pdf,
                PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES.out.concordance_detection,
                PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES.out.concordance_sample,
                PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES.out.concordance_sig_corr
            )
    }
    def ch_comparative_collected = ch_comparative_figs.collect().ifEmpty( [] )

    // ── Step 10: Gather deliverables ─────────────────────────────────────────
    // Collect per-method activity matrices and key figures into one
    // mutsig_artifacts_${cohort_id}/ folder. DBS/ID and MuSiCaL/SigDyn are
    // conditional; missing categories are simply omitted.
    def ch_activity_matrices = MERGE_ACTIVITIES_SBS.out.merged_activities
        .mix(
            MERGE_ACTIVITIES_DBS.out.merged_activities,
            MERGE_ACTIVITIES_ID.out.merged_activities
        )
    if ( ch_musical_activities ) ch_activity_matrices = ch_activity_matrices.mix( ch_musical_activities )
    if ( ch_sigdyn_activities )  ch_activity_matrices = ch_activity_matrices.mix( ch_sigdyn_activities )
    def ch_activity_collected = ch_activity_matrices.collect()

    // Signature bargraph PNG + PDF for every mutation type that ran (SBS always;
    // DBS/ID only when their matrices were non-empty).
    def ch_bargraph_figs = PLOT_BARGRAPHS_SBS.out.bargraph_png
        .mix(
            PLOT_BARGRAPHS_SBS.out.bargraph_pdf,
            PLOT_BARGRAPHS_DBS.out.bargraph_png,
            PLOT_BARGRAPHS_DBS.out.bargraph_pdf,
            PLOT_BARGRAPHS_ID.out.bargraph_png,
            PLOT_BARGRAPHS_ID.out.bargraph_pdf
        )
        .collect()

    // Interactive HTML bargraphs (not part of the gather bundle) — collected for publish.
    def ch_bargraph_htmls = PLOT_BARGRAPHS_SBS.out.bargraph_html
        .mix(
            PLOT_BARGRAPHS_DBS.out.bargraph_html,
            PLOT_BARGRAPHS_ID.out.bargraph_html
        )
        .collect()

    GATHER_MUTSIG_ARTIFACTS (
        cohort_id,
        ch_activity_collected,
        ch_bargraph_figs,
        ch_comparative_collected
    )

    emit:
    mutsig_artifacts  = GATHER_MUTSIG_ARTIFACTS.out.bundle       // tuple(cohort_id, dir)
    activity_matrices = ch_activity_collected                    // value ch: SPA(SBS/DBS/ID) + MuSiCaL + SigDyn
    musicatk_extra    = ch_musicatk_pub                          // musicatk activities/signatures (empty if off)
    bargraph_figs     = ch_bargraph_figs                         // value ch: PNG + PDF (SBS/DBS/ID)
    bargraph_htmls    = ch_bargraph_htmls                        // value ch: interactive HTML
    comparative_figs  = ch_comparative_collected                 // value ch (may be [])
    coverage_tsvs     = ch_cov_tsv_list                          // value ch: per-sample coverage TSVs
}
