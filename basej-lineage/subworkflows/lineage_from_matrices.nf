nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// LINEAGE_FROM_MATRICES  (named subworkflow, basej style)
//
// Matrices-driven lineage / phylogeny subworkflow. Instead of extracting and
// assembling NR/NV from VCFs, it consumes PRE-BUILT matrices directly and runs
// SEQUOIA phylogeny + variant placement on them:
//
//   * nr_matrix / nv_matrix  (REQUIRED) — paired NR/NV matrices (rows =
//         CHROM_POS_REF_ALT, cols = samples). Used as the single 'unfiltered'
//         scheme; NR/NV drive treemut branch assignment AND are passed to
//         SEQUOIA_VARIANT_PLACEMENT_* for variant placement.
//   * genotype_bin           (OPTIONAL, "Option D") — pre-built binarized genotype
//         matrix (0/0.5/1). When provided (non-empty sentinel), SEQUOIA_PHYLOGENY_*
//         skips internal VAF discretization; absent => /dev/null sentinel.
//   * per-sample annotated VCFs — only for the per-sample genotype tables consumed
//         by POSTPROCESS; NOT used to (re)build NR/NV.
//   * mandatory_variants_qc_status (OPTIONAL) — when set, failing mandatory variants
//         are removed from the POSTPROCESS heatmaps; absent => /dev/null sentinel.
//
// Pipeline:
//   GENOTYPE_TABLE_FROM_ANNOTATED_VCF (per sample, for postprocess)
//   SEQUOIA_PHYLOGENY_{SNV,INDEL,BOTH}        (genotype_bin + NR + NV)
//   SEQUOIA_VARIANT_PLACEMENT_{SNV,INDEL,BOTH} (tree + NR + NV)
//   POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_{SNV,INDEL,BOTH}
//   GATHER_LINEAGE_POSTPROCESS_ARTIFACTS
//
// Input validation, sequoia_phylogeny_mode parsing, MPBoot bootstrap floor, and
// file/sentinel resolution are performed by the orchestrator (main.nf) and passed
// in via take: inputs.
// ─────────────────────────────────────────────────────────────────────────────

include { GENOTYPE_TABLE_FROM_ANNOTATED_VCF          } from '../modules.nf'
include { SEQUOIA_PHYLOGENY_SNV                      } from '../modules.nf'
include { SEQUOIA_PHYLOGENY_INDEL                    } from '../modules.nf'
include { SEQUOIA_PHYLOGENY_BOTH                     } from '../modules.nf'
include { SEQUOIA_VARIANT_PLACEMENT_SNV              } from '../modules.nf'
include { SEQUOIA_VARIANT_PLACEMENT_INDEL            } from '../modules.nf'
include { SEQUOIA_VARIANT_PLACEMENT_BOTH             } from '../modules.nf'
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_SNV   } from '../modules.nf'
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_INDEL } from '../modules.nf'
include { POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_BOTH  } from '../modules.nf'
include { GATHER_LINEAGE_POSTPROCESS_ARTIFACTS       } from '../modules.nf'

workflow LINEAGE_FROM_MATRICES {

    take:
    ch_per_sample              // tuple(group, sample_name, vcf, tbi)
    cohort_id                  // val
    nr_file                    // path (pre-built NR matrix)
    nv_file                    // path (pre-built NV matrix)
    genotype_bin_file          // path (/dev/null sentinel when absent)
    mandatory_qc_file          // path (/dev/null sentinel when absent)
    gender                     // val
    run_phylo_snv              // bool
    run_phylo_indel            // bool
    run_phylo_both             // bool
    vaf_absent                 // val
    vaf_present                // val
    tree_mut_pval              // val
    keep_ancestral             // val
    create_multi_tree          // val
    genotype_conv_prob         // val
    min_pval_for_true_somatic  // val
    min_variant_reads_shared   // val
    min_vaf_shared             // val
    mpboot_path                // val
    mpboot_bootstrap           // val

    main:

    // ── Step 1: Per-sample genotype tables (CHROM_POS_REF_ALT + GT) ───────────
    // Consumed by POSTPROCESS downstream (fanned) and emitted for publishing.
    GENOTYPE_TABLE_FROM_ANNOTATED_VCF( ch_per_sample )

    GENOTYPE_TABLE_FROM_ANNOTATED_VCF.out.genotype_table
        .multiMap { it -> to_gt: it; to_pub: it }
        .set { ch_gt_fan }

    // ── Step 2: Phylogeny inputs — single 'unfiltered' scheme ────────────────
    // Three independent channel.of() so each SEQUOIA_PHYLOGENY process gets its
    // own copy (queue channels are consumed once). genotype_bin is the 5th element.
    ch_matrix_pairs_snv   = channel.of( tuple(cohort_id, 'unfiltered', nr_file, nv_file, genotype_bin_file) )
    ch_matrix_pairs_indel = channel.of( tuple(cohort_id, 'unfiltered', nr_file, nv_file, genotype_bin_file) )
    ch_matrix_pairs_both  = channel.of( tuple(cohort_id, 'unfiltered', nr_file, nv_file, genotype_bin_file) )

    SEQUOIA_PHYLOGENY_SNV (
        run_phylo_snv ? ch_matrix_pairs_snv : channel.empty(),
        gender,
        vaf_absent,
        vaf_present,
        create_multi_tree,
        mpboot_path,
        mpboot_bootstrap
    )

    SEQUOIA_PHYLOGENY_INDEL (
        run_phylo_indel ? ch_matrix_pairs_indel : channel.empty(),
        gender,
        vaf_absent,
        vaf_present,
        create_multi_tree,
        mpboot_path,
        mpboot_bootstrap
    )

    SEQUOIA_PHYLOGENY_BOTH (
        run_phylo_both ? ch_matrix_pairs_both : channel.empty(),
        gender,
        vaf_absent,
        vaf_present,
        create_multi_tree,
        mpboot_path,
        mpboot_bootstrap
    )

    // ── Step 3: Re-group per-type phylogeny outputs (sentinel when type off) ──
    // Each phylogeny output is consumed twice (grouping for placement + publish
    // emit); fan via multiMap so each consumer reads its own dedicated branch.
    if (run_phylo_snv) {
        SEQUOIA_PHYLOGENY_SNV.out.phylogeny_outputs
            .multiMap { it -> to_group: it; to_pub: it; to_gather: it }
            .set { ch_phy_snv }
        ch_snv_phy_grouped = ch_phy_snv.to_group.groupTuple(by: 0)
        ch_phy_snv_pub     = ch_phy_snv.to_pub
        ch_phy_snv_gather  = ch_phy_snv.to_gather
    } else {
        ch_snv_phy_grouped = channel.of(tuple(cohort_id, [file('/dev/null')]))
        ch_phy_snv_pub     = channel.empty()
        ch_phy_snv_gather  = channel.empty()
    }

    if (run_phylo_indel) {
        SEQUOIA_PHYLOGENY_INDEL.out.phylogeny_outputs
            .multiMap { it -> to_group: it; to_pub: it; to_gather: it }
            .set { ch_phy_indel }
        ch_indel_phy_grouped = ch_phy_indel.to_group.groupTuple(by: 0)
        ch_phy_indel_pub     = ch_phy_indel.to_pub
        ch_phy_indel_gather  = ch_phy_indel.to_gather
    } else {
        ch_indel_phy_grouped = channel.of(tuple(cohort_id, [file('/dev/null')]))
        ch_phy_indel_pub     = channel.empty()
        ch_phy_indel_gather  = channel.empty()
    }

    if (run_phylo_both) {
        SEQUOIA_PHYLOGENY_BOTH.out.phylogeny_outputs
            .multiMap { it -> to_group: it; to_pub: it; to_gather: it }
            .set { ch_phy_both }
        ch_both_phy_grouped = ch_phy_both.to_group.groupTuple(by: 0)
        ch_phy_both_pub     = ch_phy_both.to_pub
        ch_phy_both_gather  = ch_phy_both.to_gather
    } else {
        ch_both_phy_grouped = channel.of(tuple(cohort_id, [file('/dev/null')]))
        ch_phy_both_pub     = channel.empty()
        ch_phy_both_gather  = channel.empty()
    }

    // ── Step 4: NR/NV for variant placement + postprocess (independent streams) ──
    // The provided matrices ARE the 'unfiltered' matrices. multiMap fans them into
    // one copy per consumer so each gets its own (queue channels are consumed once).
    // NOTE: there are SIX consumers — the three SEQUOIA_VARIANT_PLACEMENT_* and the
    // three POSTPROCESS_* — so each needs a dedicated branch. Reusing a single
    // 'postprocess' branch across all three postprocess steps starves whichever type
    // loses the race (the single item is delivered once), so that postprocess never
    // submits and the downstream GATHER join can never complete.
    channel.of( tuple(cohort_id, nr_file, nv_file) )
        .multiMap { group, nr, nv ->
            snv:               tuple(group, nr, nv)
            indel:             tuple(group, nr, nv)
            both:              tuple(group, nr, nv)
            postprocess_snv:   tuple(group, nr, nv)
            postprocess_indel: tuple(group, nr, nv)
            postprocess_both:  tuple(group, nr, nv)
        }
        .set { ch_unfiltered }

    // Placement uses the only available scheme tree.
    def effective_placement_scheme = 'unfiltered'

    // Helper: locate the branch-length tree inside a phylogeny output subdir.
    def findPileupTree = { outputs, subdir, mutKind ->
        def out_list = outputs instanceof List ? outputs : [outputs]
        def pileup_dir = out_list.find { it.isDirectory() && it.name == subdir }
        if (!pileup_dir) return file('/dev/null')
        def names = pileup_dir.list()
        if (!names) return file('/dev/null')
        def branch_tree = null
        if (mutKind == 'both') {
            ['_both_tree_with_branch_length_selectedscheme.tree', '_snv_tree_with_branch_length_selectedscheme.tree', '_indel_tree_with_branch_length_selectedscheme.tree'].each { suf ->
                if (!branch_tree) {
                    branch_tree = names.find { it.endsWith(suf) }
                }
            }
        } else {
            branch_tree = names.find { it.endsWith("_${mutKind}_tree_with_branch_length_selectedscheme.tree") }
        }
        if (branch_tree) return pileup_dir.resolve(branch_tree)
        def tf_name = names.find { it.endsWith('.treefile') }
        return tf_name ? pileup_dir.resolve(tf_name) : file('/dev/null')
    }

    def findPhyloPlacedVariantsAll = { outputs, subdir ->
        def out_list = outputs instanceof List ? outputs : [outputs]
        def d = out_list.find { it.isDirectory() && it.name == subdir }
        if (!d) return file('/dev/null')
        def names = d.list()
        if (!names) return file('/dev/null')
        def tsv = names.find { it.endsWith('_placed_variants_all.tsv') }
        return tsv ? d.resolve(tsv) : file('/dev/null')
    }

    ch_snv_placement_input = ch_snv_phy_grouped
        .map { group, outputs ->
            def subdir = "output_snv_${effective_placement_scheme}"
            tuple(group, findPileupTree(outputs, subdir, 'snv'), findPhyloPlacedVariantsAll(outputs, subdir))
        }
        .combine(ch_unfiltered.snv, by: 0)

    ch_indel_placement_input = ch_indel_phy_grouped
        .map { group, outputs ->
            def subdir = "output_indel_${effective_placement_scheme}"
            tuple(group, findPileupTree(outputs, subdir, 'indel'), findPhyloPlacedVariantsAll(outputs, subdir))
        }
        .combine(ch_unfiltered.indel, by: 0)

    ch_both_placement_input = ch_both_phy_grouped
        .map { group, outputs ->
            def subdir = "output_both_${effective_placement_scheme}"
            tuple(group, findPileupTree(outputs, subdir, 'both'), findPhyloPlacedVariantsAll(outputs, subdir))
        }
        .combine(ch_unfiltered.both, by: 0)

    SEQUOIA_VARIANT_PLACEMENT_SNV (
        ch_snv_placement_input,
        gender,
        vaf_absent,
        vaf_present,
        tree_mut_pval,
        keep_ancestral,
        create_multi_tree,
        genotype_conv_prob,
        min_pval_for_true_somatic,
        min_variant_reads_shared,
        min_vaf_shared
    )

    SEQUOIA_VARIANT_PLACEMENT_INDEL (
        ch_indel_placement_input,
        gender,
        vaf_absent,
        vaf_present,
        tree_mut_pval,
        keep_ancestral,
        create_multi_tree,
        genotype_conv_prob,
        min_pval_for_true_somatic,
        min_variant_reads_shared,
        min_vaf_shared
    )

    SEQUOIA_VARIANT_PLACEMENT_BOTH (
        ch_both_placement_input,
        gender,
        vaf_absent,
        vaf_present,
        tree_mut_pval,
        keep_ancestral,
        create_multi_tree,
        genotype_conv_prob,
        min_pval_for_true_somatic,
        min_variant_reads_shared,
        min_vaf_shared
    )

    // Each placement output is consumed THREE times downstream: the POSTPROCESS_*
    // combine (Step 5), the GATHER placement join (Step 6), and the publish emit.
    // Referencing a process `.out` channel directly in multiple places relies on
    // broadcast-replay timing: under -resume on Nextflow >= 26 the eagerly-replayed
    // cached items are delivered to only the FIRST subscriber, so later consumers
    // see nothing and their downstream steps never submit. Fan each placement output
    // once via multiMap into dedicated branches — each consumed exactly once.
    SEQUOIA_VARIANT_PLACEMENT_SNV.out.placement_outputs
        .multiMap { it -> to_postprocess: it; to_gather: it; to_pub: it }
        .set { ch_place_snv }
    SEQUOIA_VARIANT_PLACEMENT_INDEL.out.placement_outputs
        .multiMap { it -> to_postprocess: it; to_gather: it; to_pub: it }
        .set { ch_place_indel }
    SEQUOIA_VARIANT_PLACEMENT_BOTH.out.placement_outputs
        .multiMap { it -> to_postprocess: it; to_gather: it; to_pub: it }
        .set { ch_place_both }

    // ── Step 5: Postprocess — VAF + digital heatmaps on placed variants ──────
    // Per-group genotype tables, fanned into one copy per postprocess type. As with
    // the NR/NV channel above, this is a single-item queue channel; reusing it across
    // the three POSTPROCESS_* combines would starve whichever type loses the race.
    ch_gt_fan.to_gt
        .map  { group, _sample_name, tsv -> tuple(group, tsv) }
        .groupTuple(by: 0)
        .multiMap { group, tsvs ->
            snv:   tuple(group, tsvs)
            indel: tuple(group, tsvs)
            both:  tuple(group, tsvs)
        }
        .set { ch_gt }

    ch_postprocess_snv_in = ch_place_snv.to_postprocess
        .combine(ch_unfiltered.postprocess_snv, by: 0)
        .combine(ch_gt.snv, by: 0)
        .map { group, dirs, nr, nv, gt -> tuple(group, dirs, nr, nv, gt, mandatory_qc_file, genotype_bin_file) }
    ch_postprocess_indel_in = ch_place_indel.to_postprocess
        .combine(ch_unfiltered.postprocess_indel, by: 0)
        .combine(ch_gt.indel, by: 0)
        .map { group, dirs, nr, nv, gt -> tuple(group, dirs, nr, nv, gt, mandatory_qc_file, genotype_bin_file) }
    ch_postprocess_both_in = ch_place_both.to_postprocess
        .combine(ch_unfiltered.postprocess_both, by: 0)
        .combine(ch_gt.both, by: 0)
        .map { group, dirs, nr, nv, gt -> tuple(group, dirs, nr, nv, gt, mandatory_qc_file, genotype_bin_file) }

    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_SNV   ( ch_postprocess_snv_in )
    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_INDEL ( ch_postprocess_indel_in )
    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_BOTH  ( ch_postprocess_both_in )

    // Postprocess outputs are consumed twice (GATHER join + publish emit); fan them.
    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_SNV.out.postprocess_snv
        .multiMap { it -> to_gather: it; to_pub: it }
        .set { ch_pp_snv }
    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_INDEL.out.postprocess_indel
        .multiMap { it -> to_gather: it; to_pub: it }
        .set { ch_pp_indel }
    POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_BOTH.out.postprocess_both
        .multiMap { it -> to_gather: it; to_pub: it }
        .set { ch_pp_both }

    // ── Step 6: Gather deliverables — collect every postprocess figure (PDF + rendered
    // PNG) together with the placed-variant tables and trees those figures were built
    // from, into a single lineage_artifacts_${group}/ folder. POSTPROCESS_* emit their
    // figures as INDIVIDUAL files; flatten() merges the three per-type emissions into
    // one flat file list; the three placement dirs are grouped per group and joined.
    ch_gather_postprocess = ch_pp_snv.to_gather
        .join(ch_pp_indel.to_gather, by: 0)
        .join(ch_pp_both.to_gather, by: 0)
        .map { group, snv_files, indel_files, both_files -> tuple(group, [snv_files, indel_files, both_files].flatten()) }

    ch_gather_placement = ch_place_snv.to_gather
        .join(ch_place_indel.to_gather, by: 0)
        .join(ch_place_both.to_gather, by: 0)
        // SEQUOIA_VARIANT_PLACEMENT_* emit path("output_<type>_placement_*"), and on success
        // that glob matches BOTH the output_<type>_placement_${group} directory AND the
        // output_<type>_placement_no_results sentinel file — so each type's emission is itself
        // a 2-element list. flatten() collapses it to one flat path list (path() input rejects
        // nested lists).
        .map { group, snv_dir, indel_dir, both_dir -> tuple(group, [snv_dir, indel_dir, both_dir].flatten()) }

    // Phylogeny output dirs grouped per group (so their tree PDFs are rendered
    // into the bundle). At least one variant type always runs, so this emits per group.
    ch_phylo_for_gather = ch_phy_snv_gather
        .mix(ch_phy_indel_gather, ch_phy_both_gather)
        .groupTuple(by: 0)
        .map { group, dirs -> tuple(group, dirs.flatten()) }

    ch_gather_input = ch_gather_postprocess
        .join(ch_gather_placement, by: 0)
        .join(ch_phylo_for_gather, by: 0)
        .map { group, pp_dirs, plac_dirs, phylo_dirs -> tuple(group, pp_dirs, plac_dirs, phylo_dirs) }

    GATHER_LINEAGE_POSTPROCESS_ARTIFACTS ( ch_gather_input )

    emit:
    genotype_table    = ch_gt_fan.to_pub
    phylogeny_snv     = ch_phy_snv_pub
    phylogeny_indel   = ch_phy_indel_pub
    phylogeny_both    = ch_phy_both_pub
    placement_snv     = ch_place_snv.to_pub
    placement_indel   = ch_place_indel.to_pub
    placement_both    = ch_place_both.to_pub
    heatmaps_snv      = ch_pp_snv.to_pub
    heatmaps_indel    = ch_pp_indel.to_pub
    heatmaps_both     = ch_pp_both.to_pub
    lineage_artifacts = GATHER_LINEAGE_POSTPROCESS_ARTIFACTS.out.bundle
}
