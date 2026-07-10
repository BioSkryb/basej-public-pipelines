nextflow.enable.dsl=2

// ============================================================================
// BASEJ-SOMATIC LOCAL MODULES
// ============================================================================
// Local copies of shared modules with all publishDir / publish_dir / enable_publish
// directives removed.  Output management is handled exclusively via the
// Seqera Platform  output {}  block in main.nf.
// ============================================================================

// ---------------------------------------------------------------------------
// VCF pre-processing
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Bulk germline helpers
// ---------------------------------------------------------------------------
process BULK_GET_VARIANTS_TO_FILTER {
    tag "BULK_GET_VARIANTS_TO_FILTER"

    input:
    path(input_vcf)
    path(reference)
    val(model_vcf)

    output:
    path("variants_present_bulk.txt")

    script:
    """
    if [ "${model_vcf}" = "deepvariant" ]; then
        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf} \
            | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa \
            | bcftools norm --threads ${task.cpus} -d exact \
            | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz
    else
        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf} \
            | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa \
            | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz
    fi

    bcftools index --threads ${task.cpus} -t temp.vcf.gz
    bcftools query --print-header -f '%CHROM\\_%POS\\_%REF\\_%ALT' temp.vcf.gz \
        | sort -u > variants_present_bulk.txt
    """
}

process CREATE_EMPTY_BULK_VARIANTS {
    tag "empty_bulk"

    input:
    val(_trigger)

    output:
    path("empty_bulk_variants.txt"), emit: bulk_variants

    script:
    """
    touch empty_bulk_variants.txt
    """
}

process FILTER_CHOSEN_VARIANTS_BY_BULK {
    tag "${group}"

    input:
    tuple val(group), path(chosen_variants), path(bulk_variants), path(priority_variants)

    output:
    tuple val(group), path("chosen_variants_filtered_bulk_${group}.txt"), emit: chosen_variants
    tuple val(group), path("bulk_filter_provenance_${group}.tsv"),        emit: bulk_filter_provenance

    script:
    // priority_variants (vep_priority_variants_*.tsv / merged priority list): variants in this
    // file are rescued even if they appear in bulk (status = "Rescued").
    """
    echo -e "Variant\\tRemainingAfterBulk" > bulk_filter_provenance_${group}.tsv

    # Extract rescue IDs from priority_variants TSV (col 1 = VariantId, skip header)
    if [ -s "${priority_variants}" ]; then
        awk 'NR>1{print \$1}' ${priority_variants} > rescue_ids.txt
    else
        touch rescue_ids.txt
    fi

    if [ ! -s "${bulk_variants}" ]; then
        cp ${chosen_variants} chosen_variants_filtered_bulk_${group}.txt
        awk '{print \$0"\\tPass"}' ${chosen_variants} >> bulk_filter_provenance_${group}.tsv
    else
        awk -v out="chosen_variants_filtered_bulk_${group}.txt" \\
            -v prov="bulk_filter_provenance_${group}.tsv" \\
            'FILENAME==ARGV[1] { rescue[\$0]=1; next }
             FILENAME==ARGV[2] { if (FNR==1 && (\$0 ~ /^#/ || \$0 ~ /^CHROM/)) next; bulk[\$0]=1; next }
             {
               if (\$0 in bulk) {
                 status = (\$0 in rescue) ? "Rescued" : "Fail"
               } else {
                 status = "Pass"
               }
               if (status != "Fail") print \$0 >> out;
               print \$0"\\t"status >> prov;
             }' rescue_ids.txt ${bulk_variants} ${chosen_variants}
    fi
    """
}

// ---------------------------------------------------------------------------
// Variant extraction helpers
// ---------------------------------------------------------------------------
process GET_VARIANTS_FROM_MERGED_VCF {
    tag "${group}"

    input:
    tuple val(group), path(merged_vcf)

    output:
    tuple val(group), path("all_variants_${group}.txt"), emit: all_variants

    script:
    def vcf = merged_vcf instanceof List ? merged_vcf.find { it.name.endsWith('.vcf.gz') } ?: merged_vcf[0] : merged_vcf
    """
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} | \
    awk -v OFS="_" 'BEGIN{FS="\\t"} NF>=4 { n=split(\$4,a,","); for(i=1;i<=n;i++) print \$1,\$2,\$3,a[i] }' | \
    sort -u > all_variants_${group}.txt
    """
}

process GET_LIST_POS_FROM_CHOSEN_VARIANTS {
    tag "${group}"

    input:
    tuple val(group), path(chosen_variants)

    output:
    tuple val(group), path("list_pos_variant_${group}.txt"), emit: list_pos

    script:
    """
    awk -F'_' 'NF>=2 {print \$1"\\t"\$2}' ${chosen_variants} | sort -u -k1,1 -k2,2n \
        > list_pos_variant_${group}.txt
    """
}

process FILTER_DF_NV_BY_CHOSEN_VARIANTS {
    tag "${group}"

    input:
    tuple val(group), path(df_nv), path(chosen_variants)

    output:
    tuple val(group), path("df_nv_filtered_${group}.tsv"), emit: df_nv_filtered

    script:
    """
    head -n1 ${df_nv} > df_nv_filtered_${group}.tsv
    awk -v FS="\\t" 'NR==FNR { a[\$0]; next } FNR>1 && (\$1 in a) { print \$0 }' \
        ${chosen_variants} ${df_nv} >> df_nv_filtered_${group}.tsv
    """
}

// ---------------------------------------------------------------------------
// Pileup
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// NV/NR tables
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Binomial / betabinomial filter (first pass)
// ---------------------------------------------------------------------------
process SEQUOIA_BINOM_BETABINOM_TAB_NV_NR {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(tab_nv), path(tab_nr)
    val(par_min_cov)
    val(par_max_cov)
    val(gender)

    output:
    tuple val(group), val(chr), path("*_filtering_all.txt"), emit: df_filter

    script:
    """
    wc -l ${tab_nv}
    wc -l ${tab_nr}

    head -n1 ${tab_nv} > header_nv.tsv
    head -n1 ${tab_nr} > header_nr.tsv

    cat ${tab_nv} | tail -n +2 | grep -v "NON_REF" | cat header_nv.tsv - > mat_nv.tsv
    cat ${tab_nr} | tail -n +2 | grep -v "NON_REF" | cat header_nr.tsv - > mat_nr.tsv

    wc -l mat_nv.tsv
    wc -l mat_nr.tsv

    echo -e "Launching filtering pipeline ..."
    Rscript /usr/local/bin/sequoia_filter_variants_nvnr.R \
        -v mat_nv.tsv -r mat_nr.tsv -n 1 \
        --min_cov ${par_min_cov} --max_cov ${par_max_cov} --gender_passed ${gender}

    ls *.txt | while read file; do mv \${file} filtered_${group}_${chr}_\${file}; done
    echo -e "Done"
    """
}

// ---------------------------------------------------------------------------
// Heuristic sample-level pileup filter
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Group-level table assembly
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Sequoia phylogeny (second-pass)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// VEP annotation (optional, run only when vep_cache_dir is set)
// ---------------------------------------------------------------------------
process SUBSET_MERGED_VCF_CHOSEN_VARIANTS {
    tag "${group}"

    input:
    tuple val(group), path(merged_vcf), path(chosen_variants)

    output:
    tuple val(group), path("subset_merged_${group}.vcf.gz"), path("subset_merged_${group}.vcf.gz.tbi"), emit: subset_vcf

    script:
    def vcf = merged_vcf[0]
    """
    echo -e "Creating regions from chosen_variants ..."
    awk -F'_' 'NF>=2 {print \$1"\\t"\$2-1"\\t"\$2}' ${chosen_variants} \
        | sort -k1,1 -k2,2n | uniq > regions.bed

    echo -e "Extracting variant lines in regions ..."
    bcftools view --threads ${task.cpus} -R regions.bed -H -O v ${vcf} > variants_body.vcf

    echo -e "Filtering to exact CHROM_POS_REF_ALT match ..."
    awk -v FS="\\t" -v OFS="\\t" \
        'NR==FNR { a[\$0]; next } { mid=\$1"_"\$2"_"\$4"_"\$5; if(mid in a) print \$0 }' \
        ${chosen_variants} variants_body.vcf > chosen_body.vcf

    echo -e "Building subset VCF ..."
    bcftools view --threads ${task.cpus} -h ${vcf} \
        | cat - chosen_body.vcf \
        | bcftools view --threads ${task.cpus} -Oz -o subset_merged_${group}.vcf.gz

    echo -e "Indexing subset VCF ..."
    bcftools index -t subset_merged_${group}.vcf.gz
    """
}

process SPLIT_SUBSET_VCF_BY_CHR {
    tag "${group}_${chr}"

    input:
    tuple val(group), path(subset_vcf), path(subset_vcf_tbi), val(chr)

    output:
    tuple val(group), val(chr), path("subset_${group}_${chr}.vcf.gz"), path("subset_${group}_${chr}.vcf.gz.tbi"), emit: subset_vcf_chr

    script:
    def vcf = subset_vcf[0]
    """
    bcftools view --threads ${task.cpus} -r ${chr} -Oz -o subset_${group}_${chr}.vcf.gz ${vcf}
    bcftools index -t subset_${group}_${chr}.vcf.gz
    """
}

process VEP_ANNOTATE {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(vcf_chr), path(vcf_chr_tbi)
    path(reference)
    val(species)
    val(assembly)
    path(cache_dir)

    output:
    tuple val(group), val(chr), path("vep_${group}_${chr}.vcf.gz"),          emit: vep_vcf_chr
    tuple val(group), val(chr), path("vep_${group}_${chr}_summary.html"),     emit: vep_stats

    script:
    def vcf_in    = vcf_chr[0]
    def cache_root = cache_dir.toString()
    """
    vep --input_file ${vcf_in} \
        --output_file vep_${group}_${chr}.vcf.gz \
        --stats_file vep_${group}_${chr}_summary.html \
        --format vcf \
        --vcf \
        --chr ${chr} \
        --check_existing \
        --compress_output bgzip \
        --fork ${task.cpus} \
        --species ${species} \
        --assembly ${assembly} \
        --fasta ${reference}/genome.fa \
        --cache --dir_cache ${cache_root} \
        --offline \
        --af \
        --max_af \
        --no_check_variants_order
    """
}

process SORT_INDEX_VEP_VCF {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(vep_vcf_gz)

    output:
    tuple val(group), val(chr), path("vep_sorted_${group}_${chr}.vcf.gz"), path("vep_sorted_${group}_${chr}.vcf.gz.tbi"), emit: vep_vcf_chr

    script:
    def vcf_in = vep_vcf_gz[0]
    """
    bcftools sort -Oz -o vep_sorted_${group}_${chr}.vcf.gz ${vcf_in}
    bcftools index -t vep_sorted_${group}_${chr}.vcf.gz
    """
}

process MERGE_VEP_VCF_BY_GROUP {
    tag "${group}"

    input:
    tuple val(group), path(vep_vcf_list)

    output:
    tuple val(group), path("group_${group}_vep.vcf.gz*"), emit: vep_vcf

    script:
    """
    ls vep_sorted_${group}_*.vcf.gz 2>/dev/null | sort -V > vep_list.txt
    bcftools concat -f vep_list.txt -n -Oz -o group_${group}_vep.vcf.gz
    bcftools index -t group_${group}_vep.vcf.gz
    """
}

process FILTER_VEP_GERMLINE {
    tag "${group}"

    input:
    tuple val(group), path(vep_vcf), path(vep_vcf_tbi)
    val(max_af)
    val(filter_by_existing_variation)

    output:
    tuple val(group), path("vep_sorted_${group}.vcf.gz"), path("vep_sorted_${group}.vcf.gz.tbi"), emit: sorted_vep_vcf
    tuple val(group), path("chosen_variants_postgermlinefilter_${group}.txt"),                      emit: chosen_variants
    tuple val(group), path("vep_filter_provenance_${group}.tsv"),                                   emit: filter_provenance

    script:
    def vcf      = vep_vcf[0]
    def filter_ev = filter_by_existing_variation ? 1 : 0
    """
    echo -e "Sorting VEP VCF before split-vep ..."
    bcftools sort -Oz -o vep_sorted_${group}.vcf.gz ${vcf}
    bcftools index -t vep_sorted_${group}.vcf.gz

    echo -e "Extracting VEP fields ..."
    bcftools +split-vep vep_sorted_${group}.vcf.gz \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%Existing_variation\\t%AF\\t%MAX_AF\\n' \
        -d -s worst > split_vep_raw.tsv

    echo -e "Filtering VEP VCF: filter_by_existing_variation=${filter_by_existing_variation}; exclude AF > ${max_af} ..."

    echo -e "CHROM\\tPOS\\tREF\\tALT\\tExisting_variation\\tAF\\tMAX_AF\\tFILTER_STATUS\\tFILTER_REASON" \
        > vep_filter_provenance_${group}.tsv

    # Pre-create empty chosen_variants file to ensure it always exists
    touch chosen_variants_postgermlinefilter_${group}.txt

    awk -v max_af=${max_af} -v filter_ev=${filter_ev} \
        -v prov="vep_filter_provenance_${group}.tsv" \
        -v chosen="chosen_variants_postgermlinefilter_${group}.txt" \
        'BEGIN{FS="\\t"; OFS="\\t"}
      NF>=7 {
        keep_ev=1; keep_af=1; reason="";
        if (filter_ev) {
          ev=\$5; keep_ev=0;
          if (ev=="" || ev==".") keep_ev=1;
          else {
            gsub(/&/, ",", ev); n=split(ev,a,",");
            for(i=1;i<=n;i++) if (a[i]!="" && a[i]!~/^rs[0-9]+\$/) { keep_ev=1; break }
          }
        }
        if (!((\$6=="" || \$6=="." || \$6+0<=max_af) && (\$7=="" || \$7=="." || \$7+0<=max_af))) keep_af=0;
        if (!keep_ev) reason="existing_variation_filter";
        if (!keep_af) reason=(reason=="" ? "AF_filter" : reason";AF_filter");
        status=(keep_ev && keep_af) ? "PASS" : "REMOVED";
        if (status=="PASS") reason=".";
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,status,reason >> prov;
        if (status=="PASS") {
          n=split(\$4,a,","); for(i=1;i<=n;i++) print \$1"_"\$2"_"\$3"_"a[i] >> chosen
        }
      }' split_vep_raw.tsv

    # Sort unique variants if file has content (size > 0)
    if [[ -s chosen_variants_postgermlinefilter_${group}.txt ]]; then
        sort -u chosen_variants_postgermlinefilter_${group}.txt \
            -o chosen_variants_postgermlinefilter_${group}.txt
    fi

    echo -e "Wrote \$(wc -l < chosen_variants_postgermlinefilter_${group}.txt) variants."
    echo -e "Wrote \$(wc -l < vep_filter_provenance_${group}.tsv) lines (incl. header) to provenance."
    """
}

// ---------------------------------------------------------------------------
// Filter provenance: master variant × filter-status table + diagnostic plots
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Per-sample VCF splitting and annotation
// ---------------------------------------------------------------------------
process LIST_SAMPLES_FROM_GROUP_VCF {
    tag "${group}"

    input:
    tuple val(group), path(vcf), path(tbi)

    output:
    tuple val(group), path("sample_list.txt"), emit: sample_list

    script:
    """
    bcftools query -l ${vcf} | sort > sample_list.txt
    echo "Samples found: \$(wc -l < sample_list.txt)"
    """
}

process ANNOTATE_SAMPLE_VCF {
    tag "${group}:${sample_name}"

    input:
    tuple val(group), val(sample_name),
          path(group_vcf), path(tbi),
          path(annot_table),
          path(pileup_focal)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_somatic_annotated.vcf.gz"),
          path("${sample_name}_somatic_annotated.vcf.gz.tbi"),
          emit: annotated_vcf

    script:
    """
    set -euo pipefail

    # ── A: Join both annotation sources into a single combined TSV ───────────────
    awk -F'\\t' -v smp="${sample_name}" '
    FILENAME == ARGV[1] && FNR == 1 {
        n_pileup = NF - 6
        for (i = 7; i <= NF; i++) pileup_hdr[i-6] = "SMPL_PILEUP_" \$i
        next
    }
    FILENAME == ARGV[1] && \$1 == smp && \$6 != "REF" {
        vid = \$2
        val = \$7
        for (i = 8; i <= NF; i++) val = val "\\t" \$i
        pileup_alt[vid] = val
        next
    }
    FILENAME == ARGV[1] && \$1 == smp && \$6 == "REF" {
        vid = \$2
        val = "0\\t0\\t0\\t0"
        for (i = 11; i <= NF; i++) val = val "\\t" \$i
        pileup_ref[vid] = val
        next
    }
    FILENAME == ARGV[2] && FNR == 1 {
        printf "CHROM\\tPOS\\tREF\\tALT"
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        for (i = 1; i <= n_pileup; i++) printf "\\t%s", pileup_hdr[i]
        print ""
        na_str = "NA"
        for (i = 2; i <= n_pileup; i++) na_str = na_str "\\tNA"
        next
    }
    FILENAME == ARGV[2] {
        n = split(\$1, a, "_")
        alt = a[n]; ref = a[n-1]; pos = a[n-2]
        chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
        printf "%s\\t%s\\t%s\\t%s", chrom, pos, ref, alt
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        if      (\$1 in pileup_alt) printf "\\t%s", pileup_alt[\$1]
        else if (\$1 in pileup_ref) printf "\\t%s", pileup_ref[\$1]
        else                         printf "\\t%s", na_str
        print ""
    }
    ' ${pileup_focal} ${annot_table} \
    | sort -k1,1 -k2,2n \
    | bgzip -c > combined.tsv.gz
    tabix -s1 -b2 -e2 -S1 combined.tsv.gz

    # ── B: Extract per-sample VCF with exact CHROM_POS_REF_ALT matching ──────────
    zcat combined.tsv.gz | awk 'NR>1 {print \$1"\\t"(\$2-1)"\\t"\$2}' \
        | sort -k1,1 -k2,2n | uniq > regions.bed

    zcat combined.tsv.gz | awk 'NR>1 {print \$1"_"\$2"_"\$3"_"\$4}' > focal_variants.txt

    bcftools view -s ${sample_name} -R regions.bed -Ov ${group_vcf} \
        | bcftools +setGT -- -t ./. -n 0/0 \
        > variants_tmp.vcf
    grep "^#"  variants_tmp.vcf > header.vcf
    grep -v "^#" variants_tmp.vcf > body.vcf
    awk -F'\\t' -v OFS='\\t' \
        'NR==FNR { a[\$0]; next } { if ((\$1"_"\$2"_"\$4"_"\$5) in a) print }' \
        focal_variants.txt body.vcf > chosen_body.vcf
    cat header.vcf chosen_body.vcf | bcftools view -Oz -o ${sample_name}_raw.vcf.gz
    bcftools index -t ${sample_name}_raw.vcf.gz

    # ── C: Build ##INFO header lines and -c column mapping ───────────────────────
    set +o pipefail
    HEADER_FIELDS=\$(zcat combined.tsv.gz | head -1 | cut -f5- | tr '\\t' '\\n')
    set -o pipefail

    echo "\${HEADER_FIELDS}" | awk '
        /^SMPL_PILEUP_/ { print "##INFO=<ID=" \$1 ",Number=1,Type=String,Description=\\"Per-sample pileup metric: " \$1 "\\">"; next }
                        { print "##INFO=<ID=" \$1 ",Number=1,Type=String,Description=\\"Somatic filter provenance: " \$1 "\\">"}
    ' > combined_hdr.txt

    ALL_COLS=\$(echo "\${HEADER_FIELDS}" | awk '{printf ",INFO/%s", \$1}' | sed 's/^,//')

    # ── D: Single bcftools annotate pass ─────────────────────────────────────────
    bcftools annotate \
        -a combined.tsv.gz \
        -h combined_hdr.txt \
        -c CHROM,POS,REF,ALT,\${ALL_COLS} \
        -Oz -o ${sample_name}_somatic_annotated.vcf.gz \
        ${sample_name}_raw.vcf.gz
    bcftools index -t ${sample_name}_somatic_annotated.vcf.gz

    echo "Done: \$(bcftools view -H ${sample_name}_somatic_annotated.vcf.gz | wc -l) variants annotated for ${sample_name}"
    """
}

// Identify high-confidence germline variants from statistical evidence and annotate
// the merged group VCF with twelve new INFO fields in a single bcftools annotate pass:
//
//   HIGH_CONFIDENCE_GERMLINE_FROM_STATS
//       Stage 1 — annot_table: Binom_Germline_qval_log10 > -1 AND (Binom_Rho < 0.1 OR Binom_Rho = NA)
//       Stage 2 — group VCF genotypes: > 80 % of samples carry a het call (REF + any ALT, ploidy-agnostic)
//       Values: "Yes" | "No"
//
//   VEP_AF_FILTER
//       Variants in vep_filter_provenance with FILTER_REASON == "AF_filter" (exact match).
//       Values: "AF_filter" | "Not_filtered"
//
//   VEP_GERMLINE_HIGH_CONFIDENCE
//       Subset of VEP_AF_FILTER variants where >= 80 % of cells carry a het (0/1) call.
//       Values: "Yes" | "No"
//
//   RemainingAfterBulk
//       Pass/Fail status from bulk_filter_provenance (FILTER_CHOSEN_VARIANTS_BY_BULK output).
//       Variants absent from the provenance (not chosen) receive "Not_evaluated".
//       When no bulk VCF was used (empty provenance), all variants receive "Pass".
//
//   RemainingAfterBulk_HIGH_CONFIDENCE
//       Subset of RemainingAfterBulk=Fail variants where >= 80 % of cells carry a het call
//       (REF + any ALT, ploidy-agnostic).
//       Values: "Yes" | "No"
//
// Additional outputs:
//   germline_stats_prevalence_${group}.tsv — per-variant non-REF prevalence for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants
//   af_filter_prevalence_${group}.tsv      — per-variant non-REF prevalence for VEP AF_filter variants
//   gt_vep.txt                             — AF_filter variant IDs passing the 80 % het threshold
//   bulk_fail_prevalence_${group}.tsv      — per-variant non-REF prevalence for RemainingAfterBulk=Fail variants
//   gt_bulk.txt                            — Bulk Fail variant IDs passing the 80 % het threshold
process IDENTIFY_GERMLINE_FROM_STATS {
    tag "${group}"

    input:
    tuple val(group), path(group_vcf), path(group_tbi), path(annot_table), path(vep_filter_provenance), path(bulk_filter_provenance)
    val(germline_prev_pct)

    output:
    tuple val(group),
          path("${group}_hcgermline_annotated.vcf.gz"),
          path("${group}_hcgermline_annotated.vcf.gz.tbi"),
          emit: annotated_vcf

    script:
    """
    set -euo pipefail

    # ── Stage 1: statistical filter on annot_table ────────────────────────────
    awk -F'\\t' '
    function qval_passes(v) {
        if (v == "NA" || v == "NaN") return 0
        if (v == "-Inf")             return 0
        if (v == "Inf")              return 1
        return v + 0 > -1
    }
    function rho_passes(v) {
        if (v == "NA" || v == "NaN") return 1
        if (v == "Inf" || v == "-Inf") return 0
        return v + 0 < 0.1
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "Binom_Germline_qval_log10") qcol = i
            if (\$i == "Binom_Rho")                 rcol = i
        }
        next
    }
    qcol && rcol && qval_passes(\$qcol) && rho_passes(\$rcol) { print \$1 }
    ' ${annot_table} | sort -u > candidate_ids.txt

    echo "[IDENTIFY_GERMLINE_FROM_STATS] Stage 1 candidates: \$(wc -l < candidate_ids.txt)"

    # ── Stage 2: genotype prevalence filter (> 80 % het) ─────────────────────
    n_samples=\$(bcftools query -l ${group_vcf} | wc -l)
    echo "[IDENTIFY_GERMLINE_FROM_STATS] Samples in VCF: \${n_samples}"

    touch germline_confirmed_ids_${group}.txt
    if [ -s "candidate_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' candidate_ids.txt | sort -k1,1 -k2,2n | uniq > candidate_regions.bed

        bcftools query \
            -R candidate_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${group_vcf} > candidate_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR {
            vid = \$0; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            cands[chrom "_" pos "_" ref "_" alt] = vid
            next
        }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in cands)) next
            het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            if (het / n_smp * 100 > pct) print cands[key]
        }
        ' candidate_ids.txt candidate_gts.tsv | sort -u > germline_confirmed_ids_${group}.txt
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes: \$(wc -l < germline_confirmed_ids_${group}.txt)"

    # ── HIGH_CONFIDENCE_GERMLINE_FROM_STATS: prevalence table ────────────────
    # Prevalence = cells where GT is not 0/0 (any non-REF call, including het and hom-alt).
    # candidate_gts.tsv (built in Stage 2) covers all Stage 1 candidates — a superset
    # of confirmed IDs — so no extra bcftools query is needed here.
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > germline_stats_prevalence_${group}.tsv

    if [ -s "germline_confirmed_ids_${group}.txt" ]; then
        awk -v n_smp="\${n_samples}" -F'\\t' '
        NR == FNR { conf[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in conf)) next
            nonref = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
        }
        ' germline_confirmed_ids_${group}.txt candidate_gts.tsv \
            >> germline_stats_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] germline_stats_prevalence rows: \$(tail -n+2 germline_stats_prevalence_${group}.tsv | wc -l)"

    # ── Parse VEP AF_filter variants ──────────────────────────────────────────
    if [ -s "${vep_filter_provenance}" ]; then
        awk -F'\\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if (\$i == "CHROM")         cc = i
                if (\$i == "POS")           pc = i
                if (\$i == "REF")           rc = i
                if (\$i == "ALT")           ac = i
                if (\$i == "FILTER_REASON") fc = i
            }
            next
        }
        cc && fc && \$fc == "AF_filter" { print \$cc "_" \$pc "_" \$rc "_" \$ac }
        ' ${vep_filter_provenance} | sort -u > af_filter_ids.txt
    else
        touch af_filter_ids.txt
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] AF_filter variants from VEP: \$(wc -l < af_filter_ids.txt)"

    # ── AF_filter: prevalence table + VEP_GERMLINE_HIGH_CONFIDENCE (>= 80 % het) ──
    # Prevalence = cells where GT is not REF (0/0) and not missing (./.)
    # VEP_GERMLINE_HIGH_CONFIDENCE = Yes when >= 80 % of cells have 0/1 het call.
    touch gt_vep.txt
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > af_filter_prevalence_${group}.tsv

    if [ -s "af_filter_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' af_filter_ids.txt | sort -k1,1 -k2,2n | uniq > af_filter_regions.bed

        bcftools query \
            -R af_filter_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${group_vcf} > af_filter_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR { af_ids[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in af_ids)) next
            nonref = 0; het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
            if (het / n_smp * 100 >= pct) print key > "gt_vep.txt"
        }
        ' af_filter_ids.txt af_filter_gts.tsv >> af_filter_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] VEP_GERMLINE_HIGH_CONFIDENCE=Yes: \$(wc -l < gt_vep.txt)"

    # ── Parse bulk filter provenance ──────────────────────────────────────────
    # Format: Variant<TAB>RemainingAfterBulk  (header on line 1)
    # When no bulk VCF was used the file is empty → all variants receive "Pass".
    # Variants present in the VCF but absent from the provenance were never chosen
    # variants and receive "Not_evaluated".
    if [ -s "${bulk_filter_provenance}" ]; then
        awk -F'\\t' 'NR > 1 && NF == 2 { bulk_status[\$1] = \$2 }
        END { for (v in bulk_status) print v "\\t" bulk_status[v] }' \
            ${bulk_filter_provenance} > bulk_lookup.tsv
        bulk_mode="provenance"
    else
        touch bulk_lookup.tsv
        bulk_mode="no_bulk"
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] Bulk mode: \${bulk_mode} (\$(wc -l < bulk_lookup.tsv) entries)"

    # ── Bulk Fail: prevalence table + RemainingAfterBulk_HIGH_CONFIDENCE (>= 80 % het) ──
    # Prevalence = cells where GT is not REF (0/0) and not missing (./.)
    # RemainingAfterBulk_HIGH_CONFIDENCE = Yes when >= 80 % of cells have a het call.
    touch gt_bulk.txt
    printf 'Variant\\tPrevalence_count\\tPrevalence_proportion\\n' \
        > bulk_fail_prevalence_${group}.tsv

    if [ -s "bulk_lookup.tsv" ]; then
        awk -F'\\t' '\$2 == "Fail" { print \$1 }' bulk_lookup.tsv | sort -u > bulk_fail_ids.txt
    else
        touch bulk_fail_ids.txt
    fi
    echo "[IDENTIFY_GERMLINE_FROM_STATS] RemainingAfterBulk=Fail variants: \$(wc -l < bulk_fail_ids.txt)"

    if [ -s "bulk_fail_ids.txt" ]; then
        awk -F'_' '{
            n = NF; pos = \$(n-2)
            chrom = \$1; for (j = 2; j <= n-3; j++) chrom = chrom "_" \$j
            print chrom "\\t" (pos - 1) "\\t" pos
        }' bulk_fail_ids.txt | sort -k1,1 -k2,2n | uniq > bulk_fail_regions.bed

        bcftools query \
            -R bulk_fail_regions.bed \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
            ${group_vcf} > bulk_fail_gts.tsv

        awk -v n_smp="\${n_samples}" -v pct=${germline_prev_pct} -F'\\t' '
        NR == FNR { bulk_ids[\$0] = 1; next }
        {
            key = \$1 "_" \$2 "_" \$3 "_" \$4
            if (!(key in bulk_ids)) next
            nonref = 0; het = 0
            for (i = 5; i <= NF; i++) {
                gt = \$i; gsub(/[|]/, "/", gt)
                if (gt != "0/0" && gt != "./." && gt != ".") nonref++
                n_al = split(gt, al, "/")
                has_ref = 0; has_alt = 0
                for (a = 1; a <= n_al; a++) {
                    if (al[a] == "0") has_ref = 1
                    else if (al[a] != ".") has_alt = 1
                }
                if (has_ref && has_alt) het++
            }
            printf "%s\\t%d\\t%.4f\\n", key, nonref, nonref / n_smp
            if (het / n_smp * 100 >= pct) print key > "gt_bulk.txt"
        }
        ' bulk_fail_ids.txt bulk_fail_gts.tsv >> bulk_fail_prevalence_${group}.tsv
    fi

    echo "[IDENTIFY_GERMLINE_FROM_STATS] RemainingAfterBulk_HIGH_CONFIDENCE=Yes: \$(wc -l < gt_bulk.txt)"

    # ── Build combined annotation TSV (all VCF variants, five fields) ────────
    # Columns: CHROM POS REF ALT HIGH_CONFIDENCE_GERMLINE_FROM_STATS VEP_AF_FILTER VEP_GERMLINE_HIGH_CONFIDENCE RemainingAfterBulk RemainingAfterBulk_HIGH_CONFIDENCE
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${group_vcf} | \
    awk -F'\\t' -v bulk_mode="\${bulk_mode}" '
    BEGIN {
        while ((getline line < "germline_confirmed_ids_${group}.txt") > 0) {
            vid = line; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            confirmed[chrom "_" pos "_" ref "_" alt] = 1
        }
        while ((getline line < "af_filter_ids.txt") > 0)
            af_ids[line] = 1
        while ((getline line < "gt_vep.txt") > 0)
            gt_vep_ids[line] = 1
        while ((getline line < "gt_bulk.txt") > 0)
            gt_bulk_ids[line] = 1
        while ((getline line < "bulk_lookup.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 2) bulk[a[1]] = a[2]
        }
        while ((getline line < "candidate_ids.txt") > 0) {
            vid = line; n = split(vid, a, "_")
            alt = a[n]; ref = a[n-1]; pos = a[n-2]
            chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
            stage1_ids[chrom "_" pos "_" ref "_" alt] = 1
        }
        while ((getline line < "germline_stats_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { gstat_c[a[1]] = a[2]; gstat_p[a[1]] = a[3] }
        }
        while ((getline line < "af_filter_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { vep_c[a[1]] = a[2]; vep_p[a[1]] = a[3] }
        }
        while ((getline line < "bulk_fail_prevalence_${group}.tsv") > 0) {
            n = split(line, a, "\\t")
            if (n == 3 && a[1] != "Variant") { bkf_c[a[1]] = a[2]; bkf_p[a[1]] = a[3] }
        }
    }
    {
        key = \$1 "_" \$2 "_" \$3 "_" \$4
        hcg     = (key in confirmed)   ? "Yes"       : "No"
        af      = (key in af_ids)      ? "AF_filter"  : "Not_filtered"
        vep_hc  = (key in gt_vep_ids)  ? "Yes"        : "No"
        if (bulk_mode == "no_bulk")      rab = "Pass"
        else if (key in bulk)            rab = bulk[key]
        else                             rab = "Not_evaluated"
        rab_hc  = (key in gt_bulk_ids) ? "Yes"        : "No"
        gfs     = (key in stage1_ids)  ? "Yes"        : "No"
        gs_cnt  = (key in gstat_c)     ? gstat_c[key] : "."
        gs_pct  = (key in gstat_c)     ? gstat_p[key] : "."
        vaf_cnt = (key in vep_c)       ? vep_c[key]   : "."
        vaf_pct = (key in vep_c)       ? vep_p[key]   : "."
        bk_cnt  = (key in bkf_c)       ? bkf_c[key]   : "."
        bk_pct  = (key in bkf_c)       ? bkf_p[key]   : "."
        print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" hcg "\\t" af "\\t" vep_hc "\\t" rab "\\t" rab_hc "\\t" gfs "\\t" gs_cnt "\\t" gs_pct "\\t" vaf_cnt "\\t" vaf_pct "\\t" bk_cnt "\\t" bk_pct
    }' | sort -k1,1 -k2,2n | bgzip -c > annot_combined.tsv.gz
    tabix -s1 -b2 -e2 annot_combined.tsv.gz

    # ── Write INFO header lines ───────────────────────────────────────────────
    printf '##INFO=<ID=HIGH_CONFIDENCE_GERMLINE_FROM_STATS,Number=1,Type=String,Description="Germline: Binom_Germline_qval_log10>-1 AND Binom_Rho<0.1 AND >%d%% samples het (0/1)">\n' \
        "${germline_prev_pct}" > combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER,Number=1,Type=String,Description="AF_filter if variant removed by VEP germline filter (FILTER_REASON==AF_filter exact); Not_filtered otherwise">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_GERMLINE_HIGH_CONFIDENCE,Number=1,Type=String,Description="Yes if VEP AF_filter variant with >=%d%% of cells carrying het (0/1) genotype">\n' \
        "${germline_prev_pct}" >> combined_hdr.txt
    printf '##INFO=<ID=RemainingAfterBulk,Number=1,Type=String,Description="Bulk filter status: Pass/Fail from bulk_filter_provenance; Not_evaluated if variant was not a chosen variant; Pass when no bulk VCF was used">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=RemainingAfterBulk_HIGH_CONFIDENCE,Number=1,Type=String,Description="Yes if RemainingAfterBulk=Fail variant with >=%d%% of cells carrying a het call (REF + any ALT, ploidy-agnostic)">\n' \
        "${germline_prev_pct}" >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_FROM_STATS,Number=1,Type=String,Description="Stage 1 germline: Binom_Germline_qval_log10>-1 AND (Binom_Rho<0.1 OR Binom_Rho=NA); no genotype prevalence filter">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_STATS_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=GERMLINE_STATS_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for VEP_AF_FILTER=AF_filter variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=VEP_AF_FILTER_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for VEP_AF_FILTER=AF_filter variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=BULK_FAIL_Prevalence_count,Number=1,Type=Integer,Description="Cells with non-REF genotype for RemainingAfterBulk=Fail variants; missing for all others">\n' \
        >> combined_hdr.txt
    printf '##INFO=<ID=BULK_FAIL_Prevalence_proportion,Number=1,Type=Float,Description="Proportion of cells with non-REF genotype for RemainingAfterBulk=Fail variants; missing for all others">\n' \
        >> combined_hdr.txt

    # ── Single bcftools annotate pass (twelve INFO fields) ───────────────────
    bcftools annotate \
        -a annot_combined.tsv.gz \
        -h combined_hdr.txt \
        -c CHROM,POS,REF,ALT,INFO/HIGH_CONFIDENCE_GERMLINE_FROM_STATS,INFO/VEP_AF_FILTER,INFO/VEP_GERMLINE_HIGH_CONFIDENCE,INFO/RemainingAfterBulk,INFO/RemainingAfterBulk_HIGH_CONFIDENCE,INFO/GERMLINE_FROM_STATS,INFO/GERMLINE_STATS_Prevalence_count,INFO/GERMLINE_STATS_Prevalence_proportion,INFO/VEP_AF_FILTER_Prevalence_count,INFO/VEP_AF_FILTER_Prevalence_proportion,INFO/BULK_FAIL_Prevalence_count,INFO/BULK_FAIL_Prevalence_proportion \
        -Oz -o ${group}_hcgermline_annotated.vcf.gz \
        ${group_vcf}
    bcftools index -t ${group}_hcgermline_annotated.vcf.gz

    echo "[IDENTIFY_GERMLINE_FROM_STATS] Done — final variant counts:"
    bcftools view -H ${group}_hcgermline_annotated.vcf.gz \
        | awk '
          \$8 ~ /HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes/ { hcg++ }
          \$8 ~ /GERMLINE_FROM_STATS=Yes/                 { gfs++ }
          \$8 ~ /VEP_AF_FILTER=AF_filter/                 { af++ }
          \$8 ~ /VEP_GERMLINE_HIGH_CONFIDENCE=Yes/        { vhc++ }
          \$8 ~ /RemainingAfterBulk=Pass/                 { rbp++ }
          \$8 ~ /RemainingAfterBulk=Fail/                 { rbf++ }
          \$8 ~ /RemainingAfterBulk=Not_evaluated/        { rbn++ }
          \$8 ~ /RemainingAfterBulk_HIGH_CONFIDENCE=Yes/  { rab_hc++ }
          END {
            print "  HIGH_CONFIDENCE_GERMLINE_FROM_STATS=Yes    : " hcg+0
            print "  GERMLINE_FROM_STATS=Yes                    : " gfs+0
            print "  VEP_AF_FILTER=AF_filter                    : " af+0
            print "  VEP_GERMLINE_HIGH_CONFIDENCE=Yes           : " vhc+0
            print "  RemainingAfterBulk=Pass                    : " rbp+0
            print "  RemainingAfterBulk=Fail                    : " rbf+0
            print "  RemainingAfterBulk=Not_evaluated           : " rbn+0
            print "  RemainingAfterBulk_HIGH_CONFIDENCE=Yes     : " rab_hc+0
          }'
    """
}


// Extracts a long-format prevalence TSV from the IDENTIFY_GERMLINE_FROM_STATS
// annotated VCF for three filter categories:
//
//   GERMLINE_FROM_STATS   — GT-based non-REF prevalence for all Stage-1
//                           germline candidates (Binom q-value + Rho filter).
//                           Computed from VCF genotypes; covers the FULL
//                           distribution across the prevalence threshold.
//
//   VEP_AF_filter         — non-REF prevalence (stored INFO field) for all
//                           variants flagged by the VEP AF filter.
//
//   Bulk_Fail             — non-REF prevalence (stored INFO field) for all
//                           variants that failed the bulk filter.
//
// Output TSV columns: Variant  Filter  Prevalence_proportion
//
// The TSV is consumed by PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS.

process EXTRACT_GERMLINE_PREVALENCE_TABLE {
    tag "${group}"

    input:
    tuple val(group), path(vcf), path(tbi)

    output:
    tuple val(group), path("germline_prevalence_long_${group}.tsv"), emit: prevalence_table

    script:
    """
    set -euo pipefail

    n_smp=\$(bcftools query -l ${vcf} | wc -l)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Samples: \${n_smp}"

    printf 'Variant\\tFilter\\tPrevalence_proportion\\n' \
        > germline_prevalence_long_${group}.tsv

    # ── 1. GERMLINE_FROM_STATS=Yes: GT-based prevalence (full distribution) ──
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying GERMLINE_FROM_STATS variants..."
    bcftools query \
        -i 'INFO/GERMLINE_FROM_STATS="Yes"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
        ${vcf} \
    | awk -v n_smp="\${n_smp}" -F'\\t' '{
        nonref = 0
        for (i = 5; i <= NF; i++) {
            gt = \$i; gsub(/[|]/, "/", gt)
            if (gt != "0/0" && gt != "./." && gt != ".") nonref++
        }
        printf "%s_%s_%s_%s\\tGERMLINE_FROM_STATS\\t%.6f\\n", \$1, \$2, \$3, \$4, nonref / n_smp
    }' >> germline_prevalence_long_${group}.tsv

    n_germ=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c GERMLINE || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] GERMLINE_FROM_STATS rows: \${n_germ}"

    # ── 2. VEP_AF_filter: stored prevalence field ────────────────────────────
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying VEP AF_filter variants..."
    bcftools query \
        -i 'INFO/VEP_AF_FILTER="AF_filter"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/VEP_AF_FILTER_Prevalence_proportion\\n' \
        ${vcf} \
    | awk -F'\\t' '\$5 != "." { printf "%s_%s_%s_%s\\tVEP_AF_filter\\t%s\\n", \$1, \$2, \$3, \$4, \$5 }' \
    >> germline_prevalence_long_${group}.tsv

    n_vep=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c VEP || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] VEP_AF_filter rows: \${n_vep}"

    # ── 3. Bulk_Fail: stored prevalence field ────────────────────────────────
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Querying Bulk_Fail variants..."
    bcftools query \
        -i 'INFO/RemainingAfterBulk="Fail"' \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BULK_FAIL_Prevalence_proportion\\n' \
        ${vcf} \
    | awk -F'\\t' '\$5 != "." { printf "%s_%s_%s_%s\\tBulk_Fail\\t%s\\n", \$1, \$2, \$3, \$4, \$5 }' \
    >> germline_prevalence_long_${group}.tsv

    n_bulk=\$(tail -n+2 germline_prevalence_long_${group}.tsv | grep -c Bulk || true)
    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Bulk_Fail rows: \${n_bulk}"

    echo "[EXTRACT_GERMLINE_PREVALENCE_TABLE] Total rows (excl. header): \$(tail -n+2 germline_prevalence_long_${group}.tsv | wc -l)"
    """
}


// Reads the long-format prevalence TSV produced by EXTRACT_GERMLINE_PREVALENCE_TABLE
// and generates a publication-quality histogram figure:
//
//   germline_prevalence_distributions_{group}.png
//       — one facet per filter (GERMLINE_FROM_STATS, VEP_AF_filter, Bulk_Fail)
//       — bars coloured by pass/fail relative to germline_prev_pct threshold
//       — vertical dashed line at the threshold
//       — variant counts and median annotated per panel

process PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS {
    tag "${group}"

    input:
    tuple val(group), path(prevalence_table)
    val(germline_prev_pct)

    output:
    tuple val(group), path("germline_prevalence_distributions_${group}.png"), emit: prevalence_plot

    script:
    """
    echo "${germline_prev_pct}" > threshold.txt

    Rscript /usr/local/bin/plot_germline_prevalence_distributions.R

    # Rename placeholder to group-specific filename
    mv germline_prevalence_distributions_PLACEHOLDER.png \
       germline_prevalence_distributions_${group}.png
    """
}


// One task per sample: extract three per-sample VCFs from the group VCF annotated
// by IDENTIFY_GERMLINE_FROM_STATS, each retaining a different germline call-set:
//
//   *_hcgermline_stats.vcf.gz  — HIGH_CONFIDENCE_GERMLINE_FROM_STATS = "Yes"
//                                (Binom stats filter + > 80 % het genotype prevalence)
//
//   *_hcgermline_vep.vcf.gz   — VEP_GERMLINE_HIGH_CONFIDENCE = "Yes"
//                                (VEP AF_filter variants with >= 80 % het prevalence)
//
//   *_hcgermline_bulk.vcf.gz  — RemainingAfterBulk_HIGH_CONFIDENCE = "Yes"
//                                (bulk-derived variants confirmed by >= 80 % het prevalence across cells)
process SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS {
    tag "${group}:${sample_name}"

    input:
    tuple val(group), val(sample_name), path(annotated_vcf), path(tbi)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_stats.vcf.gz"),
          path("${sample_name}_hcgermline_stats.vcf.gz.tbi"),
          emit: stats_vcf
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_vep.vcf.gz"),
          path("${sample_name}_hcgermline_vep.vcf.gz.tbi"),
          emit: vep_vcf
    tuple val(group), val(sample_name),
          path("${sample_name}_hcgermline_bulk.vcf.gz"),
          path("${sample_name}_hcgermline_bulk.vcf.gz.tbi"),
          emit: bulk_vcf

    script:
    """
    set -euo pipefail

    # ── HIGH_CONFIDENCE_GERMLINE_FROM_STATS = "Yes" ───────────────────────────
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/HIGH_CONFIDENCE_GERMLINE_FROM_STATS="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_stats.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_stats.vcf.gz

    # ── VEP_GERMLINE_HIGH_CONFIDENCE = "Yes" ─────────────────────────────────
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/VEP_GERMLINE_HIGH_CONFIDENCE="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_vep.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_vep.vcf.gz

    # ── RemainingAfterBulk_HIGH_CONFIDENCE = "Yes" (bulk-confirmed germline) ──
    bcftools view -s ${sample_name} ${annotated_vcf} \
      | bcftools filter -i 'INFO/RemainingAfterBulk_HIGH_CONFIDENCE="Yes"' \
      | bcftools view -Oz -o ${sample_name}_hcgermline_bulk.vcf.gz
    bcftools index -t ${sample_name}_hcgermline_bulk.vcf.gz

    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} stats : \$(bcftools view -H ${sample_name}_hcgermline_stats.vcf.gz | wc -l) variants"
    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} vep   : \$(bcftools view -H ${sample_name}_hcgermline_vep.vcf.gz   | wc -l) variants"
    echo "[SUBSET_MERGED_VCF_HCGERMLINE] ${sample_name} bulk  : \$(bcftools view -H ${sample_name}_hcgermline_bulk.vcf.gz  | wc -l) variants"
    """
}


// Adapted from CREATE_ADO_TABLE (modules/bioskryb/ado/create_ado_table/main.nf).
// The original module works with per-sample gVCFs and a separate baseline het-sites VCF.
// Here the input is already a single-sample germline VCF (output of
// SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS) which:
//   • contains only confirmed germline het sites → no additional het-site reference needed
//   • is a regular VCF (not gVCF) → gVCF conversion is not required
// Allele balance is computed directly from the FORMAT/AD and FORMAT/DP fields.
//
// Output df_ADO_* TSV format (no header, 10 columns):
//   CHROM  POS  REF  ALT  ALT  1  AD_ALT  DP  FREQ  PROVENANCE
// Columns 1–9 match what SUMMARIZE_ADO_INTERVALS expects.
// Column 10 (PROVENANCE) records the germline filter set ("stats", "vep", or "bulk")
// so downstream summaries can be split per filter type.

process CREATE_ADO_TABLE_FROM_GERMLINE_VCF {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(vcf), path(tbi)
    val(sample_prop)

    output:
    tuple val(sample_name), path("df_ADO_${sample_name}.tsv"), emit: ado_table

    script:
    // Provenance = the trailing suffix of sample_name added by the pipeline wiring:
    // "SAMPLE_stats", "SAMPLE_vep", or "SAMPLE_bulk" → "stats" / "vep" / "bulk".
    def provenance = sample_name.tokenize('_').last()
    """
    set -euo pipefail

    # Extract allele balance at every germline het site.
    # %AD{0} = REF depth, %AD{1} = ALT depth (first ALT; germline sites are biallelic).
    # FREQ = AD_ALT / DP  (0 when DP = 0, i.e. complete allele drop-out).
    # Column 10 records the provenance label for downstream per-set summaries.
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD{1}]\\t[%DP]\\n' \\
        ${vcf} \\
    | awk -v OFS='\\t' -v prop=${sample_prop} -v prov="${provenance}" \\
        'BEGIN { srand() }
        rand() < prop {
            ad_alt = \$5 + 0
            dp     = \$6 + 0
            freq   = (dp > 0) ? ad_alt / dp : 0
            print \$1, \$2, \$3, \$4, \$4, 1, ad_alt, dp, freq, prov
        }' \\
    > df_ADO_${sample_name}.tsv

    echo "[CREATE_ADO_TABLE_FROM_GERMLINE_VCF] ${sample_name}: \$(wc -l < df_ADO_${sample_name}.tsv) sites"
    """
}


// Consumes all per-provenance merged_ADO_*.tsv and merged_ADO_summary_*.tsv
// files produced by CONCAT_ADO_STATS / CONCAT_ADO_VEP / CONCAT_ADO_BULK and
// generates two publication-quality plots:
//
//   ADO_germline_dist.png       — allele-frequency distribution per germline
//                                 filter set: per-sample points + red mean
//                                 line/square, faceted by provenance
//   ADO_germline_summary.png    — per-sample ADO proportion in [0.2–0.8],
//                                 dashed lines connecting the same cell
//                                 across filter sets
//   ADO_germline_comparison.png — both panels combined side-by-side
//
// The R code avoids dollar-sign column access (uses [[]]) and regex end-anchors
// so the Nextflow triple-quoted string does not require excessive escaping.

process PLOT_ADO_GERMLINE_COMPARISON {
    tag "ado_germline_plots"

    input:
    path(merged_ado_tables)       // merged_ADO_{stats,vep,bulk}.tsv
    path(merged_ado_summaries)    // merged_ADO_summary_{stats,vep,bulk}.tsv

    output:
    path("ADO_germline_dist.png"),        emit: dist_plot
    path("ADO_germline_summary.png"),     emit: summary_plot
    path("ADO_germline_comparison.png"),  emit: combined_plot
    path("ADO_germline_summary.tsv"),     emit: summary_table

    script:
    """
    Rscript /usr/local/bin/plot_ado_germline_comparison.R
    """
}


// Variant of CONCAT_SUMMARY_ADO_INTERVALS that accepts a provenance label
// ("stats", "vep", or "bulk") and names all outputs with that label, so
// per-provenance summaries are kept separate in the publish directory.
//
// Inputs are the res_ADO_* TSV files collected from SUMMARIZE_ADO_INTERVALS
// after they have been split by provenance (see pipeline wiring).
//
// Output files:
//   merged_ADO_${label}.tsv        — combined allele-frequency distribution
//   ADO_plot_summary_${label}.png  — summary plot
//   merged_ADO_summary_${label}.tsv — per-sample ADO percentage

process CONCAT_SUMMARY_ADO_INTERVALS_LABELED {
    tag "${label}"

    input:
    path(df_files)
    val(label)

    output:
    path("merged_ADO_${label}.tsv"),          emit: merged_ADO
    path("ADO_plot_summary_${label}.png"),    emit: plot_ADO
    path("merged_ADO_summary_${label}.tsv"),  emit: summary_ADO

    script:
    """
    set -euo pipefail

    echo -e "File_Interval\\tFreq\\tProp" > merged_ADO_${label}.tsv

    cat res_ADO_* | cut -f1  | sed 's/\\.tsv//' > a.txt
    cat res_ADO_* | sed 's/\\[//' | sed 's/)//' | sed 's/]//' | cut -f2 > b.txt
    cat res_ADO_* | cut -f3 > c.txt
    cat res_ADO_* | cut -f4 > d.txt

    paste -d "_" a.txt b.txt > end.txt
    paste -d "\\t" end.txt c.txt d.txt >> merged_ADO_${label}.tsv

    Rscript /usr/local/bin/plot_summary_ado_intervals.R merged_ADO_${label}.tsv

    # R script writes merged_ADO_summary.tsv and ADO_plot_summary.png with hardcoded names
    mv merged_ADO_summary.tsv merged_ADO_summary_${label}.tsv
    mv ADO_plot_summary.png   ADO_plot_summary_${label}.png
    """
}

process SUMMARIZE_ADO_INTERVALS {
    tag "${sample_name}"
    
    input:
    tuple val( sample_name ), path( ado_file )
    val( cov_cutoff )
  
    output:
    path("*.png"), emit: plot
    path("res_ADO_*"), emit: df_sum

    script:
    """
    
    echo -e "Processing ${ado_file} with coverage cutoff >= ${cov_cutoff}";
    Rscript /usr/local/bin/summarize_ado_intervals.R ${ado_file} ${cov_cutoff}
    """
    
}

// ---------------------------------------------------------------------------
// NR/NV Pileup Depth Analysis (Multi-Scheme Support)
// ---------------------------------------------------------------------------

// ============================================================================
// JUNE 2026 SOMATIC UPDATE — NEW MODULES (ported from somatic_snp_indel_filtering)
// publishDir / publish-params stripped; containers set via withName in nextflow.config
// ============================================================================

// --- FILTER_MERGED_VCF_BY_ML_VERDICT
process FILTER_MERGED_VCF_BY_ML_VERDICT {
    tag "${group}"

    input:
    tuple val(group), path(merged_vcf)
    val(ml_min_pass_cells)

    output:
    tuple val(group), path("usable_merged_group_${group}.vcf.gz*"), emit: usable_vcf
    tuple val(group), path("df_nv_group_${group}.tsv"),            emit: df_nv
    tuple val(group), path("ml_verdict_stats_group_${group}.tsv.gz"), emit: ml_stats
    tuple val(group), path("ml_verdict_summary_group_${group}.txt"),  emit: ml_summary

    script:
    """

    VCF=${merged_vcf[0]}

    echo -e "Computing per-variant ML verdict counts across all cells ...";

    # Per-variant table: how many cells carry / accept / reject each variant, plus a keep/drop call.
    echo -e "variant_id\\tn_carry\\tn_pass\\tn_reject\\tfrac_reject\\tverdict" > ml_verdict_stats_group_${group}.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%MLV]\\n' "\$VCF" \\
        | awk -F'\\t' -v OFS='\\t' -v minpass=${ml_min_pass_cells} '{
              carry=0; pass=0; rej=0;
              for (i=5; i<=NF; i++) {
                  if (\$i=="0")      { carry++; pass++ }
                  else if (\$i=="1") { carry++; rej++  }
              }
              fr = (carry>0 ? rej/carry : 0);
              verdict = (pass>=minpass ? "keep" : "drop");
              print \$1"_"\$2"_"\$3"_"\$4, carry, pass, rej, fr, verdict
          }' >> ml_verdict_stats_group_${group}.tsv

    echo -e "Writing summary ...";

    total=\$(( \$(wc -l < ml_verdict_stats_group_${group}.tsv) - 1 ))
    keep=\$(awk -F'\\t' 'NR>1 && \$6=="keep"' ml_verdict_stats_group_${group}.tsv | wc -l)
    drop=\$(( total - keep ))
    pct_filtered=\$(awk -v d="\$drop" -v t="\$total" 'BEGIN{ printf "%.2f", (t>0 ? 100*d/t : 0) }')
    pct_kept=\$(awk -v k="\$keep" -v t="\$total" 'BEGIN{ printf "%.2f", (t>0 ? 100*k/t : 0) }')

    echo "==================================================================="
    echo "FILTER_MERGED_VCF_BY_ML_VERDICT  [group: ${group}]"
    echo "  ml_min_pass_cells (keep if accepted by >= this many cells): ${ml_min_pass_cells}"
    echo "  Input variants            : \$total"
    echo "  Output (usable) variants  : \$keep (\${pct_kept}% kept)"
    echo "  Filtered out variants     : \$drop (\${pct_filtered}% filtered)"
    echo "==================================================================="

    {
        printf 'group\\t%s\\n' "${group}"
        printf 'ml_min_pass_cells\\t%s\\n' "${ml_min_pass_cells}"
        printf 'total_variants\\t%s\\n' "\$total"
        printf 'usable_variants\\t%s\\n' "\$keep"
        printf 'filtered_variants\\t%s\\n' "\$drop"
        printf 'pct_kept\\t%s\\n' "\$pct_kept"
        printf 'pct_filtered\\t%s\\n' "\$pct_filtered"
    } > ml_verdict_summary_group_${group}.txt

    cat ml_verdict_summary_group_${group}.txt

    echo -e "Subsetting merged VCF to usable variants (accepted by >= ${ml_min_pass_cells} cell(s)) ...";

    bcftools view --threads ${task.cpus} -i 'N_PASS(MLV==0) >= ${ml_min_pass_cells}' "\$VCF" -Oz -o usable_merged_group_${group}.vcf.gz
    bcftools index --threads ${task.cpus} -t usable_merged_group_${group}.vcf.gz

    bgzip ml_verdict_stats_group_${group}.tsv

    echo -e "Regenerating df_nv table from the usable VCF ...";

    bcftools query -l usable_merged_group_${group}.vcf.gz > columns.txt

    cat columns.txt | tr "\\n" "\\t" | sed "s|\\t\$|\\n|" | sed "s|^|\\t|" > df_nv_group_${group}.tsv

    bcftools query --print-header -f '%CHROM\\_%POS\\_%REF\\_%ALT[\\t%AD]\\n' usable_merged_group_${group}.vcf.gz | tail -n+2 >> df_nv_group_${group}.tsv

    """
}

// --- FILTER_VEP_GERMLINE_CHR
process FILTER_VEP_GERMLINE_CHR {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(vep_vcf)
    val(max_af)
    val(filter_by_existing_variation)

    output:
    tuple val(group), val(chr), path("vep_sorted_${group}_${chr}.vcf.gz"), path("vep_sorted_${group}_${chr}.vcf.gz.tbi"), emit: sorted_vep_vcf
    tuple val(group), val(chr), path("chosen_variants_postgermlinefilter_${group}_${chr}.txt"), emit: chosen_variants
    tuple val(group), val(chr), path("vep_filter_provenance_${group}_${chr}.tsv"), emit: filter_provenance
    tuple val(group), val(chr), path("vep_variant_annotation_${group}_${chr}.tsv"), emit: variant_annotation
    tuple val(group), val(chr), path("vep_impact_verdict_table_${group}_${chr}.tsv"), emit: impact_verdict_table
    tuple val(group), val(chr), path("vep_priority_variants_${group}_${chr}.tsv"), emit: priority_variants

    script:
    def vcf = vep_vcf instanceof List ? vep_vcf[0] : vep_vcf
    def filter_ev = filter_by_existing_variation ? 1 : 0
    """
    # Sort to a temp name first, then move into place. The staged input may already be
    # named vep_sorted_${group}_${chr}.vcf.gz (it comes from SORT_INDEX_VEP_VCF), so sorting
    # directly onto that name would have bcftools open the read-only input symlink for
    # writing (Permission denied / segfault on the local executor).
    bcftools sort -Oz -o sorted_tmp_${group}_${chr}.vcf.gz ${vcf}
    mv sorted_tmp_${group}_${chr}.vcf.gz vep_sorted_${group}_${chr}.vcf.gz
    bcftools index -t vep_sorted_${group}_${chr}.vcf.gz

    echo -e "Extracting VEP fields (sharded by contig over ${task.cpus} workers) ..."
    # bcftools +split-vep has no --threads option, so parallelize across contigs with xargs -P.
    # For a single-chromosome VCF this is effectively one shard; identical result to single-threaded.
    bcftools index -s vep_sorted_${group}_${chr}.vcf.gz | cut -f1 > split_vep_contigs.txt

    cat split_vep_contigs.txt | xargs -I{} -P ${task.cpus} \\
        bcftools +split-vep vep_sorted_${group}_${chr}.vcf.gz -r {} \\
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%Existing_variation\\t%AF\\t%MAX_AF\\t%IMPACT\\t%Consequence\\t%SYMBOL\\t%CLIN_SIG\\n' \\
            -d -s worst -o split_vep_raw.shard.{}.tsv

    : > split_vep_raw.tsv
    while read -r contig; do
        cat "split_vep_raw.shard.\${contig}.tsv" >> split_vep_raw.tsv
    done < split_vep_contigs.txt
    rm -f split_vep_raw.shard.*.tsv split_vep_contigs.txt

    echo -e "Filtering VEP VCF: filter_by_existing_variation=${filter_by_existing_variation}; exclude AF, MAX_AF > ${max_af} ..."

    echo -e "CHROM\\tPOS\\tREF\\tALT\\tExisting_variation\\tAF\\tMAX_AF\\tFILTER_STATUS\\tFILTER_REASON" \\
        > vep_filter_provenance_${group}_${chr}.tsv

    echo -e "VariantId\\tImpact\\tConsequence\\tGenesAffected\\tDataBasesReported\\tVerdict" \\
        > vep_variant_annotation_${group}_${chr}.tsv

    # Ensure the chosen list exists even when no variant on this chromosome passes.
    : > chosen_variants_postgermlinefilter_${group}_${chr}.txt

    awk -v max_af=${max_af} -v filter_ev=${filter_ev} \\
        -v prov="vep_filter_provenance_${group}_${chr}.tsv" \\
        -v chosen="chosen_variants_postgermlinefilter_${group}_${chr}.txt" \\
        -v annot="vep_variant_annotation_${group}_${chr}.tsv" \\
        'BEGIN{FS="\\t"; OFS="\\t"}
      NF>=10 {
        keep_ev=1; keep_af=1; reason="";

        # --- Existing_variation filter ---
        if (filter_ev) {
          ev=\$5; keep_ev=0;
          if (ev=="" || ev==".") keep_ev=1;
          else {
            gsub(/&/, ",", ev); n=split(ev,a,",");
            for(i=1;i<=n;i++) if (a[i]!="" && a[i]!~/^rs[0-9]+\$/) { keep_ev=1; break }
          }
        }

        # --- AF / MAX_AF filter ---
        if (!((\$6=="" || \$6=="." || \$6+0<=max_af) && (\$7=="" || \$7=="." || \$7+0<=max_af))) keep_af=0;

        # --- Determine status and reason ---
        if (!keep_ev) reason="existing_variation_filter";
        if (!keep_af) reason=(reason=="" ? "AF_filter" : reason";AF_filter");
        status=(keep_ev && keep_af) ? "PASS" : "REMOVED";
        if (status=="PASS") reason=".";

        # --- Write provenance row ---
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,status,reason >> prov;

        # --- Write to chosen_variants if PASS ---
        if (status=="PASS") print \$1"_"\$2"_"\$3"_"\$4 >> chosen;

        # --- Write variant annotation row (all variants) ---
        vid=\$1"_"\$2"_"\$3"_"\$4;
        impact=(\$8=="" || \$8==".") ? "." : \$8;
        conseq=(\$9=="" || \$9==".") ? "." : \$9;
        gene=(\$10=="" || \$10==".") ? "." : \$10;
        db=(\$5=="" || \$5==".") ? "." : \$5;

        # --- Derive Verdict from CLIN_SIG (priority) or IMPACT (fallback) ---
        cs=\$11; verdict="Unknown";
        if (cs!="" && cs!=".") {
          tmp=cs; gsub(/likely_pathogenic/,"",tmp);
          if (tmp~/pathogenic/)         verdict="Pathogenic";
          else if (cs~/likely_pathogenic/) verdict="Likely_Pathogenic";
          else if (cs~/uncertain/)      verdict="VUS";
          else if (cs~/benign/)         verdict="Benign";
        } else {
          if      (impact=="HIGH")     verdict="Deleterious";
          else if (impact=="MODERATE") verdict="Possibly_Deleterious";
          else if (impact=="LOW")      verdict="Benign";
        }

        print vid,impact,conseq,gene,db,verdict >> annot;
      }' split_vep_raw.tsv

    sort -u chosen_variants_postgermlinefilter_${group}_${chr}.txt \\
        -o chosen_variants_postgermlinefilter_${group}_${chr}.txt

    echo -e "Building Impact x Verdict contingency table and priority variant subset ..."

    # Ensure ct_pairs.tsv exists even when this chromosome has no annotated variants.
    : > ct_pairs.tsv

    # Always create the priority file with its header up front, so the process output
    # exists even when this chromosome has zero passing/annotated variants.
    echo -e "VariantId\\tImpact\\tConsequence\\tGenesAffected\\tDataBasesReported\\tVerdict" \\
        > vep_priority_variants_${group}_${chr}.tsv

    # Build priority subset and collect raw impact/verdict pairs for the contingency table.
    # NOTE: keyed on FILENAME==ARGV[1] (not FNR==NR) so an EMPTY chosen file does not
    # cause every annotation row to be mis-read as a chosen entry (the FNR==NR empty-
    # first-file trap, which otherwise silently drops the priority header + all rows).
    awk -v subset="vep_priority_variants_${group}_${chr}.tsv" \\
        'BEGIN{FS="\\t"; OFS="\\t"}
      FILENAME==ARGV[1] { chosen[\$1]=1; next }
      FNR==1  { next }
      {
        vid=\$1; impact=\$2; verdict=\$6;
        print impact, verdict >> "ct_pairs.tsv"
        if ((vid in chosen) &&
            ((impact=="HIGH" && verdict!="Benign") ||
             ((verdict=="Pathogenic" || verdict=="Likely_Pathogenic") && (impact=="HIGH" || impact=="MODERATE")))) {
          print \$0 >> subset;
        }
      }' chosen_variants_postgermlinefilter_${group}_${chr}.txt vep_variant_annotation_${group}_${chr}.tsv

    # Build contingency table using shell sort (avoids gawk-only asort())
    cut -f1 ct_pairs.tsv | sort -u > impact_labels.txt
    cut -f2 ct_pairs.tsv | sort -u > verdict_labels.txt

    awk -F'\\t' \\
        'FILENAME==ARGV[1] { verd[++nv]=\$1; next }
         FILENAME==ARGV[2] { imp[++ni]=\$1;  next }
         { ct[\$1,\$2]++;   next }
         END {
           printf "Impact"
           for (j=1; j<=nv; j++) printf "\\t%s", verd[j]
           printf "\\n"
           for (i=1; i<=ni; i++) {
             printf "%s", imp[i]
             for (j=1; j<=nv; j++) {
               key = imp[i] SUBSEP verd[j]
               printf "\\t%s", (key in ct ? ct[key] : 0)
             }
             printf "\\n"
           }
         }' verdict_labels.txt impact_labels.txt ct_pairs.tsv \\
        > vep_impact_verdict_table_${group}_${chr}.tsv

    echo -e "Wrote \$(wc -l < chosen_variants_postgermlinefilter_${group}_${chr}.txt) variants to chosen_variants_postgermlinefilter_${group}_${chr}.txt"
    echo -e "Wrote \$(wc -l < vep_filter_provenance_${group}_${chr}.tsv) lines (incl. header) to vep_filter_provenance_${group}_${chr}.tsv"
    echo -e "Wrote \$(wc -l < vep_priority_variants_${group}_${chr}.tsv) lines (incl. header) to vep_priority_variants_${group}_${chr}.tsv"
    """
}

// --- FILTER_CHOSEN_VARIANTS_BY_VEP_CHR
process FILTER_CHOSEN_VARIANTS_BY_VEP_CHR {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(chosen_variants), path(vep_pass_variants)

    output:
    tuple val(group), val(chr), path("chosen_variants_filtered_vep_${group}_${chr}.txt"), emit: chosen_variants

    script:
    """
    if [ ! -s "${vep_pass_variants}" ]; then
      : > chosen_variants_filtered_vep_${group}_${chr}.txt
    else
      awk 'NR==FNR { if (FNR==1 && (\$0 ~ /^#/ || \$0 ~ /^CHROM/)) next; keep[\$0]=1; next } (\$0 in keep)' ${vep_pass_variants} ${chosen_variants} > chosen_variants_filtered_vep_${group}_${chr}.txt
    fi

    # Defensive: ensure the filtered output exists even if an awk short-circuit left it unwritten.
    [ -f chosen_variants_filtered_vep_${group}_${chr}.txt ] || : > chosen_variants_filtered_vep_${group}_${chr}.txt
    """
}

// --- CONCAT_VEP_PERCHR_OUTPUTS
process CONCAT_VEP_PERCHR_OUTPUTS {
    tag "${group}"

    input:
    tuple val(group), path(chosen_chr_files), path(provenance_chr_files), path(priority_chr_files)

    output:
    tuple val(group), path("chosen_variants_filtered_vep_${group}.txt"), emit: chosen_variants
    tuple val(group), path("vep_filter_provenance_${group}.tsv"), emit: filter_provenance
    tuple val(group), path("vep_priority_variants_${group}.tsv"), emit: priority_variants

    script:
    """
    set -euo pipefail

    # ── chosen_variants: headerless CHROM_POS_REF_ALT list, sorted + unique ──────
    : > chosen_variants_filtered_vep_${group}.txt
    cat ${chosen_chr_files} >> chosen_variants_filtered_vep_${group}.txt || true
    sort -u chosen_variants_filtered_vep_${group}.txt -o chosen_variants_filtered_vep_${group}.txt

    # ── filter_provenance: one header (from any shard) + every data row ──────────
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tExisting_variation\\tAF\\tMAX_AF\\tFILTER_STATUS\\tFILTER_REASON" \\
        > vep_filter_provenance_${group}.tsv
    for f in ${provenance_chr_files} ; do
        tail -n +2 "\$f" >> vep_filter_provenance_${group}.tsv
    done

    # ── priority_variants: one header (from any shard) + every data row ──────────
    echo -e "VariantId\\tImpact\\tConsequence\\tGenesAffected\\tDataBasesReported\\tVerdict" \\
        > vep_priority_variants_${group}.tsv
    for f in ${priority_chr_files} ; do
        tail -n +2 "\$f" >> vep_priority_variants_${group}.tsv
    done

    echo -e "Group ${group}: \$(wc -l < chosen_variants_filtered_vep_${group}.txt) chosen variants, \$(( \$(wc -l < vep_filter_provenance_${group}.tsv) - 1 )) provenance rows, \$(( \$(wc -l < vep_priority_variants_${group}.tsv) - 1 )) priority variants"
    """
}

// --- MERGE_MANDATORY_PRIORITY_VARIANTS
process MERGE_MANDATORY_PRIORITY_VARIANTS {
    tag "${group}"

    input:
    tuple val(group), path(priority_variants), path(mandatory_variants)

    output:
    tuple val(group), path("merged_priority_variants_${group}.tsv"), emit: merged_priority_variants

    script:
    """
    # Initialize output with header from priority_variants, or create one if the file is empty/sentinel
    if [ -s "${priority_variants}" ]; then
        head -n1 "${priority_variants}" > merged_priority_variants_${group}.tsv
        tail -n+2 "${priority_variants}" >> merged_priority_variants_${group}.tsv
    else
        printf 'VariantId\\tImpact\\tConsequence\\tGenesAffected\\tDataBasesReported\\tVerdict\\n' \\
            > merged_priority_variants_${group}.tsv
    fi

    # Merge mandatory variants: add those not already present in priority_variants
    if [ -s "${mandatory_variants}" ]; then
        awk -v OFS="\\t" -v outfile="merged_priority_variants_${group}.tsv" \\
            'FILENAME==ARGV[1] { if (FNR>1) seen[\$1]=1; next }
             NF>0 {
               vid=\$1;
               if (vid!="" && !(vid in seen)) {
                 print vid, ".", ".", ".", ".", "MandatoryVariant" >> outfile;
                 seen[vid]=1
               }
             }' \\
            merged_priority_variants_${group}.tsv "${mandatory_variants}"

        n_vep=\$(tail -n+2 "${priority_variants}" 2>/dev/null | wc -l || echo 0)
        n_total=\$(tail -n+2 merged_priority_variants_${group}.tsv | wc -l)
        n_added=\$(( n_total - n_vep ))
        echo "VEP priority variants: \${n_vep}"
        echo "Mandatory variants added (not already in VEP set): \${n_added}"
        echo "Total merged priority variants: \${n_total}"
    else
        n_total=\$(tail -n+2 merged_priority_variants_${group}.tsv | wc -l)
        echo "No mandatory_variants file provided; keeping VEP priority variants only: \${n_total}"
    fi
    """
}


// --- VARIANT_PREFILTER_TABLE
process VARIANT_PREFILTER_TABLE {
    tag "${group}"

    input:
    tuple val(group), path(pileup_agg_partials), path(sequoia_filt)

    output:
    tuple val(group), path("variant_prefilter_table_${group}.tsv"), emit: prefilter_table

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/variant_prefilter_table.R "${group}"
    """
}

// --- VARIANT_PREFILTER_TABLE_PERCHR
process VARIANT_PREFILTER_TABLE_PERCHR {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(pileup_table)

    output:
    tuple val(group), path("prefilter_pileup_agg_${group}_${chr}.tsv"), emit: agg

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/variant_prefilter_table_perchr.R "${group}" "${chr}"
    """
}

// --- VAF_SPLIT_VARIANTS_HEXBIN
process VAF_SPLIT_VARIANTS_HEXBIN {
    tag "${group}"

    input:
    tuple val(group), path(nv), path(nr), path(master), path(mandatory)
    val(rho_thr)
    val(log10q_thr)
    val(bin_nv)
    val(bin_nv_shared)
    val(bin_nr)
    val(bin_vaf)
    val(bin_vaf_shared_anchor)
    val(bin_vaf_singleton)
    val(first_rho_thr)
    val(first_log10q_thr)

    output:
    tuple val(group), path("${group}_*"), emit: all_outputs
    tuple val(group),
          path("${group}_binary_matrix_HQRoundStatisticalFiltered.tsv"),
          path("${group}_binary_matrix_HQRoundStatisticalFilteredPlusQCFiltered.tsv"),
          path("${group}_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv"),
          path("${group}_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv"),
          path("${group}_binary_matrix_unfiltered.tsv"),
          emit: cascade_binary
    tuple val(group), path("${group}_hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.pdf"), emit: hexbin, optional: true
    tuple val(group),
          path("${group}_NV_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv"),
          path("${group}_NR_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv"),
          path("${group}_binary_matrix_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv"),
          emit: phylogeny_matrices, optional: true
    tuple val(group), path("${group}_per_cell_SNV_INDEL_counts_by_stage.pdf"), emit: per_cell_counts, optional: true

    script:
    """
    set -euo pipefail

    # --mandatory is optional: only pass it when a non-empty priority/mandatory
    # file was staged (ch_priority_variants emits /dev/null as a size-0 sentinel
    # when neither VEP priority nor params.mandatory_variants are set).
    mand_arg=""
    if [ -s "${mandatory}" ]; then
        mand_arg="--mandatory ${mandatory}"
    fi

    Rscript /usr/local/bin/build_vaf_split_variants_hexbin.R \\
        --nv     ${nv} \\
        --nr     ${nr} \\
        --master ${master} \\
        --outdir . \\
        --prefix ${group}_ \\
        --rho-thr                ${rho_thr} \\
        --log10q-thr             ${log10q_thr} \\
        --bin-nv                 ${bin_nv} \\
        --bin-nv-shared          ${bin_nv_shared} \\
        --bin-nr                 ${bin_nr} \\
        --bin-vaf                ${bin_vaf} \\
        --bin-vaf-shared-anchor  ${bin_vaf_shared_anchor} \\
        --bin-vaf-singleton      ${bin_vaf_singleton} \\
        --first-rho-thr          ${first_rho_thr} \\
        --first-log10q-thr       ${first_log10q_thr} \\
        \$mand_arg
    """
}

// --- VARIANT_FILTER_FUNNEL
process VARIANT_FILTER_FUNNEL {
    tag "${group}"

    input:
    tuple val(group), path(master_table), path(all_variants),
          path(hqstat_binary), path(qc_binary), path(depth_binary), path(phylo_binary),
          path(unfiltered_binary)

    output:
    tuple val(group), path("variant_filter_tracking_${group}.tsv"), emit: filter_tracking
    tuple val(group), path("variant_filter_report_${group}.md"),    emit: filter_report
    tuple val(group), path("variant_filter_plot_${group}.pdf"),     emit: filter_plot
    tuple val(group), path("variant_filter_plot_${group}.png"),     emit: filter_plot_png

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/variant_filter_funnel.R "${group}"
    """
}


// --- BUILD_FOCAL_PILEUP
process BUILD_FOCAL_PILEUP {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(pileup_chr), path(vcf_annotation_table), path(sample_list)

    output:
    tuple val(group), path("pileup_focal_variants_${group}_${chr}.tsv"), emit: pileup_focal_chr

    script:
    def g = group
    def c = chr
    """
    set -euo pipefail

    # Focal (VEP-PASS) VariantIds for THIS chromosome only. The trailing "_" anchors the
    # chr token (e.g. ^chr1_ matches chr1_... but not chr10_...). || true: a chr with no
    # focal variants yields an empty list (header-only output, merged away by concat).
    tail -n +2 ${vcf_annotation_table} | cut -f1 > all_focal_ids.txt
    grep "^${c}_" all_focal_ids.txt | sort > focal_variant_ids_${g}_${c}.txt || true
    echo "Focal variants on ${c}: \$(wc -l < focal_variant_ids_${g}_${c}.txt)"

    OUT=pileup_focal_variants_${g}_${c}.tsv
    head -1 ${pileup_chr} > "\$OUT"

    # Three-lookup awk pass over this chr's pileup (see module comments for the lookup logic).
    awk -F'\\t' -v OFS='\\t' '
        NR==FNR {
            alt_ids[\$1] = 1
            n = split(\$1, a, "_")
            pos = a[n-2]
            chrom = a[1]; for (j=2; j<=n-3; j++) chrom = chrom "_" a[j]
            key = chrom "_" pos
            pos_to_vid[key] = (key in pos_to_vid) ? pos_to_vid[key] "|" \$1 : \$1
            next
        }
        FNR==1 { next }
        \$2 in alt_ids { print; next }
        (\$3 "_" \$4) in pos_to_vid && \$6 == "REF" {
            nvids = split(pos_to_vid[\$3 "_" \$4], vids, "|")
            for (vi = 1; vi <= nvids; vi++) { \$2 = vids[vi]; print }
            next
        }
        (\$3 "_" \$4) in pos_to_vid {
            \$6 = "REF"; \$7 = 0; \$8 = 0; \$9 = 0; \$10 = 0
            nvids = split(pos_to_vid[\$3 "_" \$4], vids, "|")
            for (vi = 1; vi <= nvids; vi++) { \$2 = vids[vi]; print }
        }
    ' focal_variant_ids_${g}_${c}.txt ${pileup_chr} >> "\$OUT"

    # Deduplicate per (sample x VariantId): ALT row wins; REF row only when no ALT exists.
    awk -F'\\t' -v OFS='\\t' '
        NR==1 { hdr=\$0; next }
        \$6 != "REF" { alt[\$1 SUBSEP \$2] = \$0; next }
                     { ref[\$1 SUBSEP \$2] = \$0 }
        END {
            print hdr
            for (k in alt) print alt[k]
            for (k in ref) if (!(k in alt)) print ref[k]
        }
    ' "\$OUT" > "\${OUT}.dedup" && mv "\${OUT}.dedup" "\$OUT"

    # Synthetic 0-rows for (global sample x this-chr focal-variant) combos absent from the
    # pileup. Global samples come from sample_list (== pileup sample set). Cols 7-22 = 0,
    # cols 23+ = NA, matching the pileup column schema.
    awk -F'\\t' -v OFS='\\t' '
        FILENAME==ARGV[1] && FNR==1 {
            ncols = NF
            filler = "0"
            for (i=8;  i<=22;    i++) filler = filler OFS "0"
            for (i=23; i<=ncols; i++) filler = filler OFS "NA"
            next
        }
        FILENAME==ARGV[2]            { samples[\$1]=1; next }
        FILENAME==ARGV[3]            { focal[\$1]=1;   next }
        FILENAME==ARGV[4] && FNR>1   { present[\$1 SUBSEP \$2]=1; next }
        END {
            for (smp in samples) {
                for (vid in focal) {
                    if (!((smp SUBSEP vid) in present)) {
                        n = split(vid, a, "_")
                        alt=a[n]; ref=a[n-1]; pos=a[n-2]
                        chrom=a[1]; for(j=2;j<=n-3;j++) chrom=chrom"_"a[j]
                        print smp, vid, chrom, pos, ref, alt, filler
                    }
                }
            }
        }
    ' "\$OUT" ${sample_list} focal_variant_ids_${g}_${c}.txt "\$OUT" >> "\$OUT"

    echo "${c} focal pileup rows: \$(tail -n +2 "\$OUT" | wc -l)"
    """
}

// --- CONCAT_FOCAL_PILEUP
process CONCAT_FOCAL_PILEUP {
    tag "${group}"

    input:
    tuple val(group), path(chr_files)

    output:
    tuple val(group), path("pileup_focal_variants_${group}.tsv"), emit: pileup_focal

    script:
    """
    set -euo pipefail
    files="${chr_files}"
    first=\$(echo \$files | tr ' ' '\\n' | head -1)
    head -1 "\$first" > pileup_focal_variants_${group}.tsv
    for f in \$files; do
        tail -n +2 "\$f"
    done >> pileup_focal_variants_${group}.tsv
    echo "Merged \$(echo \$files | wc -w) per-chr files; total focal rows: \$(tail -n +2 pileup_focal_variants_${group}.tsv | wc -l)"
    """
}

// --- ANNOTATE_VCF_SCHEME_MEMBERSHIP
process ANNOTATE_VCF_SCHEME_MEMBERSHIP {
    tag "${group}:${sample_name}"

    input:
    tuple val(group), val(sample_name),
          path(vcf), path(tbi),
          path(membership)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_somatic_annotated_schemes.vcf.gz"),
          path("${sample_name}_somatic_annotated_schemes.vcf.gz.tbi"),
          emit: annotated_vcf

    script:
    """
    set -euo pipefail

    # ── A: membership TSV (VariantId + 5 flags) -> CHROM POS REF ALT + flags ─────
    # VariantId is CHROM_POS_REF_ALT: ALT = last '_' field, REF = n-1, POS = n-2,
    # CHROM = everything before (handles names like chrUn_gl000220).
    awk -F'\\t' '
    NR == 1 {
        printf "#CHROM\\tPOS\\tREF\\tALT"
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        print ""; next
    }
    {
        n = split(\$1, a, "_")
        alt = a[n]; ref = a[n-1]; pos = a[n-2]
        chrom = a[1]; for (j = 2; j <= n-3; j++) chrom = chrom "_" a[j]
        printf "%s\\t%s\\t%s\\t%s", chrom, pos, ref, alt
        for (i = 2; i <= NF; i++) printf "\\t%s", \$i
        print ""
    }' ${membership} > scheme_ann_unsorted.tsv

    # header row first, then coordinate-sorted body
    head -n1 scheme_ann_unsorted.tsv > scheme_ann.tsv
    tail -n +2 scheme_ann_unsorted.tsv | sort -k1,1 -k2,2n >> scheme_ann.tsv
    bgzip -c scheme_ann.tsv > scheme_ann.tsv.gz
    tabix -s1 -b2 -e2 -S1 scheme_ann.tsv.gz

    # ── B: ##INFO header lines (one per scheme flag) ─────────────────────────────
    cat > scheme_hdr.txt <<'HDR'
##INFO=<ID=Scheme_unfiltered,Number=1,Type=Integer,Description="Variant is in the 'unfiltered' NR/NV matrix scheme (1) or not (0)">
##INFO=<ID=Scheme_HQStat,Number=1,Type=Integer,Description="Variant is in the HQRoundStatisticalFiltered scheme (1) or not (0)">
##INFO=<ID=Scheme_HQStat_QC,Number=1,Type=Integer,Description="Variant is in the HQRoundStatisticalFilteredPlusQCFiltered scheme (1) or not (0)">
##INFO=<ID=Scheme_HQStat_QC_Depth,Number=1,Type=Integer,Description="Variant is in the HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered scheme (1) or not (0)">
##INFO=<ID=Scheme_HQStat_QC_Depth_Phylo,Number=1,Type=Integer,Description="Variant is in the HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny scheme (1) or not (0)">
HDR

    # ── C: annotate ──────────────────────────────────────────────────────────────
    bcftools annotate \
        -a scheme_ann.tsv.gz \
        -h scheme_hdr.txt \
        -c CHROM,POS,REF,ALT,INFO/Scheme_unfiltered,INFO/Scheme_HQStat,INFO/Scheme_HQStat_QC,INFO/Scheme_HQStat_QC_Depth,INFO/Scheme_HQStat_QC_Depth_Phylo \
        -Oz -o ${sample_name}_somatic_annotated_schemes.vcf.gz \
        ${vcf}
    bcftools index -t ${sample_name}_somatic_annotated_schemes.vcf.gz

    echo "Done: tagged \$(bcftools view -H ${sample_name}_somatic_annotated_schemes.vcf.gz | wc -l) variants for ${sample_name}"
    """
}

// --- MERGE_ANNOTATED_SAMPLE_VCFS
process MERGE_ANNOTATED_SAMPLE_VCFS {
    tag "${group}"

    input:
    tuple val(group), path(vcfs), path(tbis)

    output:
    tuple val(group),
          path("${group}_merged_annotated_schemes.vcf.gz"),
          path("${group}_merged_annotated_schemes.vcf.gz.tbi"),
          emit: merged_vcf

    script:
    """
    set -euo pipefail

    # Deterministic input order for bcftools merge.
    ls *_somatic_annotated_schemes.vcf.gz | sort > vcf_list.txt
    n_vcf=\$(wc -l < vcf_list.txt)
    echo "Merging \${n_vcf} per-sample VCF(s) for group ${group}"

    if [ "\${n_vcf}" -eq 1 ]; then
        # A single cell: nothing to merge, just carry it through.
        cp "\$(cat vcf_list.txt)" merged_raw.vcf.gz
    else
        bcftools merge -l vcf_list.txt -Oz -o merged_raw.vcf.gz
    fi
    bcftools index -t merged_raw.vcf.gz

    # Strip per-sample SMPL_PILEUP_* INFO (cannot merge meaningfully across samples).
    rm_fields=\$(bcftools view -h merged_raw.vcf.gz \
        | grep -oE 'ID=SMPL_PILEUP_[A-Za-z0-9_]+' \
        | sed 's#^ID=#INFO/#' | sort -u | paste -sd, -)

    if [ -n "\${rm_fields}" ]; then
        echo "Stripping per-sample INFO: \${rm_fields}"
        bcftools annotate -x "\${rm_fields}" \
            -Oz -o ${group}_merged_annotated_schemes.vcf.gz merged_raw.vcf.gz
    else
        cp merged_raw.vcf.gz ${group}_merged_annotated_schemes.vcf.gz
    fi
    bcftools index -t ${group}_merged_annotated_schemes.vcf.gz

    echo "Done: \$(bcftools view -H ${group}_merged_annotated_schemes.vcf.gz | wc -l) variants x \$(bcftools query -l ${group}_merged_annotated_schemes.vcf.gz | wc -l) samples"
    """
}

// --- SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH
// Subset the merged, scheme-annotated group VCF to ONLY the variants in the
// HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered scheme (INFO/Scheme_HQStat_QC_Depth=1,
// added by ANNOTATE_VCF_SCHEME_MEMBERSHIP). All sample columns / INFO fields preserved.
process SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH {
    tag "${group}"

    input:
    tuple val(group), path(vcf), path(tbi)

    output:
    tuple val(group),
          path("${group}_HQStat_QC_Depth.vcf.gz"),
          path("${group}_HQStat_QC_Depth.vcf.gz.tbi"),
          emit: subset_vcf

    script:
    """
    set -euo pipefail

    bcftools view -i 'INFO/Scheme_HQStat_QC_Depth=1' \
        -Oz -o ${group}_HQStat_QC_Depth.vcf.gz \
        ${vcf}
    bcftools index -t ${group}_HQStat_QC_Depth.vcf.gz

    echo "Done: \$(bcftools view -H ${group}_HQStat_QC_Depth.vcf.gz | wc -l) HQStat_QC_Depth variants (of \$(bcftools view -H ${vcf} | wc -l) total) for group ${group}"
    """
}

// --- EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF
process EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF {
    tag "${group}__${sample_name}"

    input:
    tuple val(group), val(sample_name), path(vcf), path(tbi)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_nr.txt"),
          path("${sample_name}_nv.txt"),
          path("${sample_name}_gt.txt"),
          path("${sample_name}_variant_ids.txt"),
          emit: per_sample_vectors

    script:
    """
    set -euo pipefail

    # NR: locus depth denominator = total HQ fragments at the POSITION.
    # SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION is annotated per VCF record (CHROM_POS_REF_ALT) and
    # is NA->0 for alleles this sample does not carry. Depth is a property of the locus, not the
    # allele, so resolve NR by CHROM_POS: every record at a position takes the MAX position depth
    # seen across all co-located records. Two-pass awk preserves VCF record order so nr.txt stays
    # row-aligned with nv/gt/variant_ids below. (Matches group-level Mat_NR at multiallelic/indel loci.)
    bcftools query \\
        -f '%CHROM\\t%POS\\t%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_POSITION\\n' \\
        "${vcf}" \\
        | awk -F'\\t' '
            { key = \$1"_"\$2; v = (\$3+0); order[NR] = key; if (v > maxd[key]) maxd[key] = v; n = NR }
            END { for (i = 1; i <= n; i++) printf "%d\\n", maxd[order[i]] }
          ' > "${sample_name}_nr.txt"

    # NV: HQ ALT-supporting fragments (forward + reverse)
    bcftools query \\
        -f '%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_F\\t%INFO/SMPL_PILEUP_NUM_FRAGMENTS_HQ_MQ_BQ_R\\n' \\
        "${vcf}" \\
        | awk -F'\\t' '{printf "%d\\n", (\$1+0)+(\$2+0)}' > "${sample_name}_nv.txt"

    # GT: genotype call — used for per-sample presence/absence
    bcftools query \\
        -f '[%GT]\\n' \\
        "${vcf}" > "${sample_name}_gt.txt"

    # Variant IDs — CHROM_POS_REF_ALT (identical across all samples in the group)
    bcftools query \\
        -f '%CHROM\\_%POS\\_%REF\\_%ALT\\n' \\
        "${vcf}" > "${sample_name}_variant_ids.txt"
    """
}

// --- CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF
process CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF {
    tag "${group}"

    input:
    tuple val(group), val(sample_names),
          path(nr_files),
          path(nv_files),
          path(gt_files),
          path(variant_ids_files),
          path(master)
    val(rho_thr)
    val(log10q_thr)
    val(bin_nv)
    val(bin_nv_shared)
    val(bin_nr)
    val(bin_vaf)
    val(bin_vaf_shared_anchor)
    val(bin_vaf_singleton)

    output:
    tuple val(group),
          path("NR_annotated_vcf_${group}_*.tsv"),
          path("NV_annotated_vcf_${group}_*.tsv"),
          emit: nr_nv_matrices
    tuple val(group),
          path("matrix_scheme_summary_${group}.tsv"),
          emit: matrix_scheme_summary
    tuple val(group),
          path("matrix_per_sample_summary_${group}.tsv"),
          emit: matrix_per_sample_summary
    tuple val(group),
          path("scheme_membership_${group}.tsv"),
          emit: scheme_membership
    // First-round (Binom) Rho-vs-qval hexbin: build_vaf_split produces it here because this
    // module passes the FULL master (which carries Binom_Rho/Binom_Germline_qval), unlike
    // VAF_SPLIT_VARIANTS_HEXBIN whose prefilter --master lacks those columns. Always emitted
    // (empty placeholder if absent) so the MultiQC join never drops the group.
    tuple val(group),
          path("${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png"),
          emit: first_round_hexbin

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Sample list + variant IDs from the pre-extracted per-sample vectors
    # ─────────────────────────────────────────────────────────────────────────
    ls *_nr.txt | sort | sed 's/_nr\\.txt//' > sample_list.txt
    n_samples=\$(wc -l < sample_list.txt)
    echo "Samples: \${n_samples}"
    if [ "\${n_samples}" -lt 1 ]; then
        echo "CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF: sample_list.txt has no samples (group=${group})" >&2
        exit 1
    fi

    # Variant IDs — all *_variant_ids.txt files are identical; use any one
    cp \$(ls *_variant_ids.txt | head -1) variant_ids_all.txt
    echo "Total variants: \$(wc -l < variant_ids_all.txt)"

    # ─────────────────────────────────────────────────────────────────────────
    # Helper: assemble NR / NV matrices from a variant ID list
    # ─────────────────────────────────────────────────────────────────────────
    build_matrix() {
        local ids=\$1 label=\$2
        {
            printf ""
            while IFS= read -r s; do printf "\\t%s" "\$s"; done < sample_list.txt
            printf "\\n"
            paste "\$ids" \$(while IFS= read -r s; do echo "\${s}_nr.txt"; done < sample_list.txt)
        } > NR_annotated_vcf_${group}_\${label}.tsv

        {
            printf ""
            while IFS= read -r s; do printf "\\t%s" "\$s"; done < sample_list.txt
            printf "\\n"
            paste "\$ids" \$(while IFS= read -r s; do echo "\${s}_nv.txt"; done < sample_list.txt)
        } > NV_annotated_vcf_${group}_\${label}.tsv
    }

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Full (pre-gate) unfiltered NR / NV matrices. These are the cascade
    # R input AND the base for the gated `unfiltered` output (gate applied in STEP 4).
    # ─────────────────────────────────────────────────────────────────────────
    build_matrix    variant_ids_all.txt unfiltered

    paste variant_ids_all.txt \$(while IFS= read -r s; do echo "\${s}_nr.txt"; done < sample_list.txt) \\
        > all_nr_table.txt

    # Pristine full (pre-gate) copies: the cascade R runs on these.
    cp NR_annotated_vcf_${group}_unfiltered.tsv full_NR_unfiltered.tsv
    cp NV_annotated_vcf_${group}_unfiltered.tsv full_NV_unfiltered.tsv

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: VAF-split cascade schemes — same R logic as VAF_SPLIT_VARIANTS_HEXBIN.
    # Run on the full unfiltered NR/NV, then map the 4 cascade NR/NV outputs to this
    # module's NR_annotated_vcf_${group}_<label>.tsv naming and stage each scheme's
    # anchor-aware binary matrix (BIN_*) for the per-sample summary (STEP 6).
    # The R script also writes VAF/plot side-outputs under cascade/; they are
    # intentionally not captured by this module's output globs.
    # ─────────────────────────────────────────────────────────────────────────
    mkdir -p cascade
    Rscript /usr/local/bin/build_vaf_split_variants_hexbin.R \\
        --nv     full_NV_unfiltered.tsv \\
        --nr     full_NR_unfiltered.tsv \\
        --master ${master} \\
        --outdir cascade \\
        --prefix ${group}_ \\
        --rho-thr                ${rho_thr} \\
        --log10q-thr             ${log10q_thr} \\
        --bin-nv                 ${bin_nv} \\
        --bin-nv-shared          ${bin_nv_shared} \\
        --bin-nr                 ${bin_nr} \\
        --bin-vaf                ${bin_vaf} \\
        --bin-vaf-shared-anchor  ${bin_vaf_shared_anchor} \\
        --bin-vaf-singleton      ${bin_vaf_singleton}

    CASCADE_LABELS="HQRoundStatisticalFiltered HQRoundStatisticalFilteredPlusQCFiltered HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny"

    for label in \${CASCADE_LABELS}; do
        if [ ! -s "cascade/${group}_NR_\${label}.tsv" ] || [ ! -s "cascade/${group}_NV_\${label}.tsv" ]; then
            echo "CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF: cascade output missing for scheme '\${label}' (group=${group})." >&2
            exit 1
        fi
        cp "cascade/${group}_NR_\${label}.tsv" "NR_annotated_vcf_${group}_\${label}.tsv"
        cp "cascade/${group}_NV_\${label}.tsv" "NV_annotated_vcf_${group}_\${label}.tsv"
        # Anchor-aware binary matrix for the per-sample summary — reuse the one
        # build_vaf_split already wrote (identical "called cell" basis to VAF_SPLIT's per_cell plot).
        if [ ! -s "cascade/${group}_binary_matrix_\${label}.tsv" ]; then
            echo "CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF: cascade binary matrix missing for scheme '\${label}' (group=${group})." >&2
            exit 1
        fi
        cp "cascade/${group}_binary_matrix_\${label}.tsv" "BIN_annotated_vcf_${group}_\${label}.tsv"
    done

    # First-round (Binom) hexbin: build_vaf_split wrote it under cascade/ (the full master
    # carries Binom_Rho/Binom_Germline_qval). Publish it for the MultiQC report. Empty
    # placeholder if absent (e.g. master without Binom_* cols) so the output always exists.
    if [ -s "cascade/${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png" ]; then
        cp "cascade/${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png" "${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png"
        cp "cascade/${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.pdf" "${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.pdf" 2>/dev/null || true
    else
        echo "[CREATE_NR_NV_MATRICES] WARN: first-round hexbin not produced (master lacks Binom_Rho/Binom_Germline_qval?); emitting empty placeholder." >&2
        : > "${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png"
    fi

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4: Unfiltered-only gate — applied to the `unfiltered` NR/NV files only.
    # Keep rows that BOTH (a) have >= unfiltered_nr_min_sample_pct% of samples with
    # NR >= unfiltered_nr_min_depth, AND (b) carry at least one variant read (NV>0 in
    # >= 1 cell). (b) drops variants that are NV==0 in every cell — covered but with no
    # alt evidence anywhere — matching the variant universe of the group-level matrices.
    # The cascade schemes (STEP 3) already exclude all-NV-zero variants and are not gated here.
    # ─────────────────────────────────────────────────────────────────────────
    # (a) NR depth/cohort gate
    awk -F'\\t' -v pct="${params.unfiltered_nr_min_sample_pct}" -v depth="${params.unfiltered_nr_min_depth}" -v n="\${n_samples}" '
        n > 0 {
            pass_count = 0
            for (i=2; i<=NF; i++) if (\$i+0 >= depth) pass_count++
            if (pass_count / n * 100 >= pct) print \$1
        }
    ' all_nr_table.txt | sort -u > union_nr_gate.txt
    # (b) NV evidence: variants with NV>0 in >= 1 cell (drop variants that are NV==0 in every cell)
    awk -F'\\t' 'FNR>1 { s=0; for (i=2; i<=NF; i++) s += (\$i+0); if (s > 0) print \$1 }' \\
        full_NV_unfiltered.tsv | sort -u > union_nv_evidence.txt
    # unfiltered keep-set = (a) AND (b)
    comm -12 union_nr_gate.txt union_nv_evidence.txt > union_unfiltered_nr_gate.txt

    for kind in NR NV; do
        src=\${kind}_annotated_vcf_${group}_unfiltered.tsv
        cp "\$src" "_presieve_\${kind}.tsv"
        awk 'NR==FNR {ids[\$1]=1; next} FNR==1 || (\$1 in ids)' \\
            union_unfiltered_nr_gate.txt "\$src" > "\${src}.tmp"
        mv "\${src}.tmp" "\$src"
        if [ ! -s "\$src" ] && [ -s "_presieve_\${kind}.tsv" ]; then
            echo "[CREATE_NR_NV_MATRICES] WARN: \${src} empty after unfiltered NR gate; restoring header row only (no variant rows)." >&2
            head -n1 "_presieve_\${kind}.tsv" > "\$src"
        fi
        rm -f "_presieve_\${kind}.tsv"
    done

    # Binary matrix for the (gated) unfiltered scheme: subset the full unfiltered binary
    # that build_vaf_split wrote (NV>0 presence) to the gated variant set. Dropped by the
    # scheme-summary plot, but emitted for completeness / consistency with the other schemes.
    awk 'NR==FNR {ids[\$1]=1; next} FNR==1 || (\$1 in ids)' \\
        union_unfiltered_nr_gate.txt "cascade/${group}_binary_matrix_unfiltered.tsv" \\
        > "BIN_annotated_vcf_${group}_unfiltered.tsv"

    if [ ! -s NR_annotated_vcf_${group}_unfiltered.tsv ] || [ ! -s NV_annotated_vcf_${group}_unfiltered.tsv ]; then
        echo "CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF: unfiltered NR/NV matrix file(s) still empty after gate + header restore (group=${group})." >&2
        exit 1
    fi

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Cohort SNV / indel counts per scheme (from each NR matrix's VariantId column)
    # ─────────────────────────────────────────────────────────────────────────
    SCHEMES="unfiltered \${CASCADE_LABELS}"

    count_snv_indel_from_matrix() {
        local mat=\$1
        if [ ! -s "\${mat}" ]; then printf "0\\t0\\n"; return; fi
        tail -n +2 "\${mat}" | cut -f1 | awk '{
            n = split(\$0, a, "_")
            ref = a[n-1]; alt = a[n]
            if (length(ref) == 1 && length(alt) == 1) snv++
            else indel++
        } END { printf "%d\\t%d\\n", snv+0, indel+0 }'
    }

    printf "scheme\\tNumberOfSNVs\\tNumberOfIndels\\n" > matrix_scheme_summary_${group}.tsv
    for label in \${SCHEMES}; do
        counts=\$(count_snv_indel_from_matrix "NR_annotated_vcf_${group}_\${label}.tsv")
        printf "%s\\t%s\\n" "\${label}" "\${counts}" >> matrix_scheme_summary_${group}.tsv
    done
    echo "[matrix_scheme_summary] Done:"
    cat matrix_scheme_summary_${group}.tsv

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 6: Per-sample SNV / INDEL counts per scheme from the anchor-aware BINARY
    # matrices (called cell = value 1) — same basis as VAF_SPLIT_VARIANTS_HEXBIN's
    # per-cell plot, so the two figures' per-cell counts agree by construction.
    # ─────────────────────────────────────────────────────────────────────────
    printf "scheme\\tsample\\tNumberOfSNVs\\tNumberOfIndels\\n" > matrix_per_sample_summary_${group}.tsv

    for bin_file in BIN_annotated_vcf_${group}_*.tsv; do
        label=\$(basename "\${bin_file}" | sed 's/BIN_annotated_vcf_${group}_//' | sed 's/\\.tsv\$//')
        awk -F'\\t' -v scheme="\${label}" '
        NR==1 {
            for (i=2; i<=NF; i++) samples[i]=\$i
            ncols=NF; next
        }
        {
            n = split(\$1, a, "_")
            ref = a[n-1]; alt = a[n]
            is_snv = (length(ref)==1 && length(alt)==1)
            for (i=2; i<=ncols; i++) {
                if (\$i == 1) {
                    if (is_snv) snv[i]++
                    else        indel[i]++
                }
            }
        }
        END {
            for (i=2; i<=ncols; i++)
                printf "%s\\t%s\\t%d\\t%d\\n", scheme, samples[i], snv[i]+0, indel[i]+0
        }' "\${bin_file}" >> matrix_per_sample_summary_${group}.tsv
    done

    echo "[matrix_per_sample_summary] Done: \$(tail -n +2 matrix_per_sample_summary_${group}.tsv | wc -l) rows"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 7: Per-variant scheme-membership table — VariantId + a 0/1 flag per scheme
    # (short aliases). Universe = all extracted variants (variant_ids_all). A flag is 1
    # iff the variant is a row in that scheme's matrix. Consumed by
    # ANNOTATE_VCF_SCHEME_MEMBERSHIP to tag the per-sample VCFs for scheme subsetting.
    # ─────────────────────────────────────────────────────────────────────────
    tail -n +2 NR_annotated_vcf_${group}_unfiltered.tsv                                                       | cut -f1 | sort -u > _sch_unfiltered.txt
    tail -n +2 NR_annotated_vcf_${group}_HQRoundStatisticalFiltered.tsv                                       | cut -f1 | sort -u > _sch_hqstat.txt
    tail -n +2 NR_annotated_vcf_${group}_HQRoundStatisticalFilteredPlusQCFiltered.tsv                         | cut -f1 | sort -u > _sch_qc.txt
    tail -n +2 NR_annotated_vcf_${group}_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered.tsv        | cut -f1 | sort -u > _sch_depth.txt
    tail -n +2 NR_annotated_vcf_${group}_HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny.tsv | cut -f1 | sort -u > _sch_phylo.txt

    printf "VariantId\\tScheme_unfiltered\\tScheme_HQStat\\tScheme_HQStat_QC\\tScheme_HQStat_QC_Depth\\tScheme_HQStat_QC_Depth_Phylo\\n" > scheme_membership_${group}.tsv
    awk -F'\\t' '
        FILENAME=="_sch_unfiltered.txt" { u[\$1]=1; next }
        FILENAME=="_sch_hqstat.txt"     { h[\$1]=1; next }
        FILENAME=="_sch_qc.txt"         { q[\$1]=1; next }
        FILENAME=="_sch_depth.txt"      { d[\$1]=1; next }
        FILENAME=="_sch_phylo.txt"      { p[\$1]=1; next }
        {
            v=\$1
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", v, (v in u), (v in h), (v in q), (v in d), (v in p)
        }
    ' _sch_unfiltered.txt _sch_hqstat.txt _sch_qc.txt _sch_depth.txt _sch_phylo.txt variant_ids_all.txt \\
        >> scheme_membership_${group}.tsv
    echo "[scheme_membership] Done: \$(tail -n +2 scheme_membership_${group}.tsv | wc -l) variants"
    """
}


// --- COMPARE_HQSTAT_QC_REDUNDANCY
process COMPARE_HQSTAT_QC_REDUNDANCY {
    tag "${group}"

    input:
    tuple val(group), path(master_table), path(scheme_membership)

    output:
    tuple val(group), path("hqstat_vs_qc_contingency_${group}.tsv"), emit: contingency
    tuple val(group), path("hqstat_vs_qc_redundancy_${group}.md"),   emit: report
    tuple val(group), path("hqstat_vs_qc_heatmap_${group}.pdf"),     emit: heatmap_pdf
    tuple val(group), path("hqstat_vs_qc_heatmap_${group}.png"),     emit: heatmap_png

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/compare_hqstat_qc_redundancy.R "${group}" "${master_table}" "${scheme_membership}"
    cat hqstat_vs_qc_redundancy_${group}.md
    """
}

// --- QUANTIFY_QC_FILTER_INFLUENCE
process QUANTIFY_QC_FILTER_INFLUENCE {
    tag "${group}"

    input:
    tuple val(group), path(focal_pileup), path(scheme_membership)

    output:
    tuple val(group), path("qc_filter_influence_${group}.tsv"), emit: influence_tsv
    tuple val(group), path("qc_filter_influence_${group}.pdf"), emit: influence_pdf
    tuple val(group), path("qc_filter_influence_${group}.png"), emit: influence_png

    script:
    """
    set -euo pipefail

    # ── counts: variant-level leave-one-out + cell-variant fail combinations ──────
    # First file = scheme_membership (build HQStat set); second = focal pileup (restricted to HQStat).
    awk -F'\\t' '
    NR == FNR {
        if (FNR == 1) { for (i = 1; i <= NF; i++) if (\$i == "Scheme_HQStat") sh = i;
                        if (!sh) { print "ERROR: scheme_membership missing Scheme_HQStat" > "/dev/stderr"; exit 1 } ; next }
        if (\$sh == 1) hq[\$1] = 1
        next
    }
    FNR == 1 {
        for (i = 1; i <= NF; i++) {
            if (\$i == "VariantId")           cv = i
            if (\$i == "AS_Filter")           ca = i
            if (\$i == "PropClipped_Filter")  cc = i
            if (\$i == "BPPos_Filter")        cb = i
        }
        if (!cv || !ca || !cc || !cb) { print "ERROR: focal pileup missing a required column" > "/dev/stderr"; exit 1 }
        next
    }
    {
        vid = \$cv
        if (!(vid in hq)) next   # only variants past the 2nd-pass statistical filter (HQStat)
        asP = (\$ca == "Pass"); clP = (\$cc == "Pass"); bpP = (\$cb == "Pass")

        # variant survives artifact filtering if >=1 cell passes all three (AS & PropClipped & BPPos)
        if (asP && clP && bpP) full[vid] = 1
        if (clP && bpP)        dAS[vid]  = 1   # drop AS
        if (asP && bpP)        dCL[vid]  = 1   # drop PropClipped
        if (asP && clP)        dBP[vid]  = 1   # drop BPPos

        # cell-variant fail breakdown among the 3 artifact filters (evaluable cells only)
        if (\$ca != "NA") {
            ncv++
            if (!asP) gAS++
            if (!clP) gCL++
            if (!bpP) gBP++
            k = ((!asP) ? "AS" : "") ((!clP) ? (asP ? "Clip" : ",Clip") : "") ((!bpP) ? ((asP && clP) ? "BPPos" : ",BPPos") : "")
            if (k == "") k = "pass_all_3"
            combo[k]++
        }
    }
    END {
        nf = na = ncl = nbp = 0
        for (v in full) nf++; for (v in dAS) na++; for (v in dCL) ncl++; for (v in dBP) nbp++

        f = "qc_filter_influence_${group}.tsv"
        printf "section\\tkey\\tvalue\\n"                              > f
        printf "variant_loo\\tcandidate_full\\t%d\\n",       nf        >> f
        printf "variant_loo\\tAS_removes\\t%d\\n",           na - nf   >> f
        printf "variant_loo\\tPropClipped_removes\\t%d\\n",  ncl - nf  >> f
        printf "variant_loo\\tBPPos_removes\\t%d\\n",        nbp - nf  >> f
        printf "cellvariant_gross\\tAS_fail\\t%d\\n",          gAS+0   >> f
        printf "cellvariant_gross\\tPropClipped_fail\\t%d\\n", gCL+0   >> f
        printf "cellvariant_gross\\tBPPos_fail\\t%d\\n",       gBP+0   >> f
        printf "cellvariant_gross\\tevaluable_cellvariants\\t%d\\n", ncv+0 >> f
        for (k in combo) printf "cellvariant_combo\\t%s\\t%d\\n", k, combo[k] >> f
    }
    ' ${scheme_membership} ${focal_pileup}

    echo "[quantify_qc_filter_influence] counts:"; cat qc_filter_influence_${group}.tsv

    # ── figure (R written at runtime; no baked script) ───────────────────────────
    cat > plot_qc_influence.R <<'RPLOT'
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(egg) })
if (file.exists("/usr/local/bin/theme_ohchibi_pubr.R")) source("/usr/local/bin/theme_ohchibi_pubr.R") else theme_ohchibi_pubr <- ggplot2::theme_bw
d <- fread("qc_filter_influence_${group}.tsv", sep = "\t", header = TRUE, data.table = FALSE)

# Panel A — variant-level marginal removal per filter
loo <- d[d[, 1] == "variant_loo" & d[, 2] != "candidate_full", ]
loo[, 2] <- sub("_removes", "", loo[, 2])
loo[, 2] <- factor(loo[, 2], levels = c("AS", "PropClipped", "BPPos"))
pA <- ggplot(loo, aes(reorder(get("key"), get("value")), get("value"))) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = get("value")), hjust = -0.2, fontface = "bold", size = 4) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Variants removed by each artifact filter (leave-one-out)",
       subtitle = "${group}  -  # variants rescued if that filter alone were dropped",
       x = NULL, y = "variants removed") +
  theme_ohchibi_pubr()

# Panel B — cell-variant fail combinations among the 3 artifact filters
cmb <- d[d[, 1] == "cellvariant_combo", c(2, 3)]
names(cmb) <- c("combo", "n"); cmb[, 2] <- as.numeric(cmb[, 2])
cmb <- cmb[order(cmb[, 2]), ]
cmb[, 1] <- factor(cmb[, 1], levels = cmb[, 1])
pB <- ggplot(cmb, aes(get("combo"), get("n"))) +
  geom_col(fill = "#7B3F00", width = 0.7) +
  geom_text(aes(label = format(get("n"), big.mark = ",")), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Cell-variant filter-fail combinations (AS / PropClipped / BPPos)",
       subtitle = "how many evaluable cell-variants fail each set of filters", x = NULL, y = "cell-variants") +
  theme_ohchibi_pubr()

g <- egg::ggarrange(pA, pB, ncol = 1, heights = c(1, 1.4), draw = FALSE)
ggsave("qc_filter_influence_${group}.png", g, width = 9, height = 8, dpi = 200, bg = "white")
ggsave("qc_filter_influence_${group}.pdf", g, width = 9, height = 8)
RPLOT
    Rscript plot_qc_influence.R
    echo "[quantify_qc_filter_influence] Done."
    """
}

// --- MANDATORY_VARIANTS_QC_STATUS
process MANDATORY_VARIANTS_QC_STATUS {
    tag "${group}"

    input:
    tuple val(group), path(master_table), path(scheme_membership), path(mandatory)

    output:
    tuple val(group), path("mandatory_variants_qc_status_${group}.tsv"), emit: status_table

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/mandatory_variants_qc_status.R "${group}" "${master_table}" "${scheme_membership}" "${mandatory}"
    cat mandatory_variants_qc_status_${group}.tsv
    """
}

// --- COLLECT_DOWNSTREAM_ARTIFACTS
process COLLECT_DOWNSTREAM_ARTIFACTS {
    tag "${group}"

    input:
    tuple val(group), path(vaf_split_outputs), path(master_table), path(mandatory_status), path(scheme_membership)

    output:
    tuple val(group), path("downstream_inputs_${group}"), emit: bundle

    script:
    """
    set -euo pipefail
    dest="downstream_inputs_${group}"
    mkdir -p "\$dest"

    # VAF-split cascade matrices (Binary / VAF / NV / NR, all schemes). Globs are tolerant:
    # any pattern with no match is skipped rather than aborting.
    cp ${group}_NR_*.tsv             "\$dest"/ 2>/dev/null || true
    cp ${group}_NV_*.tsv             "\$dest"/ 2>/dev/null || true
    cp ${group}_VAF_*.tsv            "\$dest"/ 2>/dev/null || true
    cp ${group}_binary_matrix*.tsv   "\$dest"/ 2>/dev/null || true

    # Master table + mandatory-variants QC verdict + per-variant scheme-membership flags
    cp "${master_table}"      "\$dest"/
    cp "${mandatory_status}"  "\$dest"/
    cp "${scheme_membership}" "\$dest"/

    n=\$(ls -1 "\$dest" | wc -l)
    echo "[collect_downstream_artifacts] ${group}: collected \${n} files into \$dest"
    ls -1 "\$dest"
    if [ "\${n}" -lt 3 ]; then
        echo "[collect_downstream_artifacts] WARN: only \${n} files collected (expected matrices + master + mandatory)" >&2
    fi
    """
}

// --- EXPLORE_QC_FILTER_THRESHOLDS
process EXPLORE_QC_FILTER_THRESHOLDS {
    tag "${group}"

    input:
    tuple val(group), path(pileup_files), path(scheme_membership)
    val(cutoff_as)
    val(cutoff_prop_clipped)
    val(cutoff_prop_bp_under)
    val(cutoff_prop_bp_upper)
    val(cutoff_sd_indiv)
    val(cutoff_mad_indiv)
    val(cutoff_sd_both)
    val(cutoff_mad_both)
    val(cutoff_sd_extreme)
    val(cutoff_mad_extreme)

    output:
    tuple val(group),
          path("qc_threshold_distributions_${group}.pdf"),
          path("qc_threshold_*_${group}.png"),
          path("qc_threshold_summary_${group}.tsv"),
          emit: threshold_report

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/explore_qc_filter_thresholds.R \\
        "${group}" "${scheme_membership}" \\
        ${cutoff_as} ${cutoff_prop_clipped} ${cutoff_prop_bp_under} ${cutoff_prop_bp_upper} \\
        ${cutoff_sd_indiv} ${cutoff_mad_indiv} ${cutoff_sd_both} ${cutoff_mad_both} \\
        ${cutoff_sd_extreme} ${cutoff_mad_extreme}
    """
}

// --- COHORT_METRICS_SUMMARY
process COHORT_METRICS_SUMMARY {
    tag "${group}"

    input:
    tuple val(group), path(funnel_tracking), path(scheme_summary), path(contingency), path(qc_influence)

    output:
    tuple val(group), path("cohort_metrics_${group}.tsv"),  emit: metrics_tsv
    tuple val(group), path("cohort_metrics_${group}.json"), emit: metrics_json

    script:
    """
    set -euo pipefail

    flat="cohort_metrics_${group}.tsv"
    printf "section\\tkey\\tvalue\\n"        > "\$flat"
    printf "meta\\tgroup\\t${group}\\n"     >> "\$flat"

    # ── funnel: per-stage counts ─────────────────────────────────────────────────
    if [ -s "${funnel_tracking}" ]; then
        awk -F'\\t' 'NR>1 {
            printf "funnel\\t%s_n_input\\t%s\\n",    \$1, \$2
            printf "funnel\\t%s_n_passing\\t%s\\n",  \$1, \$3
            printf "funnel\\t%s_n_filtered\\t%s\\n", \$1, \$4
        }' "${funnel_tracking}" >> "\$flat"
    fi

    # ── scheme sizes: SNV / indel / total per scheme ─────────────────────────────
    if [ -s "${scheme_summary}" ]; then
        awk -F'\\t' 'NR>1 {
            printf "scheme\\t%s_snv\\t%s\\n",   \$1, \$2
            printf "scheme\\t%s_indel\\t%s\\n", \$1, \$3
            printf "scheme\\t%s_total\\t%d\\n", \$1, (\$2 + \$3)
        }' "${scheme_summary}" >> "\$flat"
    fi

    # ── HQStat-vs-QC redundancy metrics ──────────────────────────────────────────
    if [ -s "${contingency}" ]; then
        awk -F'\\t' 'NR>1 { printf "redundancy\\t%s\\t%s\\n", \$1, \$2 }' "${contingency}" >> "\$flat"
    fi

    # ── QC artifact-filter influence (already section/key/value) ─────────────────
    if [ -s "${qc_influence}" ]; then
        awk -F'\\t' 'NR>1 { printf "qc_%s\\t%s\\t%s\\n", \$1, \$2, \$3 }' "${qc_influence}" >> "\$flat"
    fi

    # ── nested JSON from the flat table ──────────────────────────────────────────
    awk -F'\\t' '
    NR==1 { next }
    {
        sec=\$1; k=\$2; v=\$3
        if (!(sec in seen)) { order[++ns]=sec; seen[sec]=1 }
        cnt[sec]++; K[sec,cnt[sec]]=k; V[sec,cnt[sec]]=v
    }
    END {
        print "{"
        for (i=1; i<=ns; i++) {
            sec=order[i]
            printf "  \\"%s\\": {\\n", sec
            for (j=1; j<=cnt[sec]; j++) {
                val=V[sec,j]
                sep=(j<cnt[sec]) ? "," : ""
                if (val ~ /^-?[0-9]+(\\.[0-9]+)?\$/) printf "    \\"%s\\": %s%s\\n", K[sec,j], val, sep
                else                                 printf "    \\"%s\\": \\"%s\\"%s\\n", K[sec,j], val, sep
            }
            printf "  }%s\\n", (i<ns) ? "," : ""
        }
        print "}"
    }' "\$flat" > "cohort_metrics_${group}.json"

    n_metrics=\$(( \$(wc -l < "\$flat") - 1 ))
    echo "[cohort_metrics_summary] Wrote \${n_metrics} metrics for ${group}"
    """
}

// --- ANALYZE_NR_NV_PILEUP_DEPTH
process ANALYZE_NR_NV_PILEUP_DEPTH {
    tag "${group}"

    // NR/NV paths must be `path` inputs so Nextflow stages them into the task dir (S3/work-dir safe).
    // Do not pass matrix paths only inside a string `val`; R would see non-local paths and fail to read.
    input:
    tuple val(group), path(nr_matrix_files), path(nv_matrix_files), val(scheme_labels_csv)
    val(min_sample_pct)
    val(depth_thresholds_csv)

    output:
    tuple val(group), path("pileup_depth_cohort_counts_${group}.tsv"),     emit: cohort_counts
    tuple val(group), path("pileup_depth_per_sc_${group}.tsv"),            emit: per_sc_counts
    tuple val(group), path("pileup_depth_per_sc_median_${group}.tsv"),     emit: per_sc_median
    tuple val(group), path("pileup_sparsity_vaf_${group}.tsv"),            emit: sparsity_tsv
    tuple val(group), path("pileup_sparsity_validation_${group}.tsv"),    emit: sparsity_validation
    tuple val(group), path("pileup_sparsity_hexbin_${group}.pdf"),         emit: sparsity_pdf
    tuple val(group), path("pileup_sparsity_hexbin_${group}_page*.png"),   emit: sparsity_png

    script:
    def nrs = nr_matrix_files instanceof List ? nr_matrix_files : [nr_matrix_files]
    def nvs = nv_matrix_files instanceof List ? nv_matrix_files : [nv_matrix_files]
    def labels = scheme_labels_csv.split(',').collect { it.trim() }.findAll { it }
    if (nrs.size() != nvs.size() || nrs.size() != labels.size()) {
        error "ANALYZE_NR_NV_PILEUP_DEPTH: mismatched nr (${nrs.size()}), nv (${nvs.size()}), labels (${labels.size()}) for group=${group}"
    }
    def argv_tail = (0..<nrs.size()).collect { i -> "\"${nrs[i]}\" \"${nvs[i]}\" \"${labels[i]}\"" }.join(' ')
    """
    set -euo pipefail
    Rscript /usr/local/bin/analyze_nr_nv_pileup_depth.R \\
        "${group}" "${min_sample_pct}" "${depth_thresholds_csv}" ${argv_tail}

    # PNG copies of the (multi-page) sparsity hexbin PDF — one PNG per page.
    pdftoppm -png -r 150 "pileup_sparsity_hexbin_${group}.pdf" "pileup_sparsity_hexbin_${group}_page"
    """
}

// --- PLOT_MATRIX_SCHEME_SUMMARY
process PLOT_MATRIX_SCHEME_SUMMARY {
    tag "${group}"

    input:
    tuple val(group), path(scheme_tsv), path(per_sample_tsv), path(upstream_per_sample_tsv)

    output:
    tuple val(group), path("matrix_scheme_summary_${group}.pdf"),      emit: scheme_summary_pdf
    tuple val(group), path("matrix_scheme_summary_${group}_page*.png"), emit: scheme_summary_png

    script:
    """

    set -euo pipefail
    Rscript /usr/local/bin/plot_matrix_scheme_summary.R \
        "${group}" "${scheme_tsv}" "${per_sample_tsv}" "${upstream_per_sample_tsv}"

    # PNG copies of the (multi-page) summary PDF — one PNG per page.
    pdftoppm -png -r 150 "matrix_scheme_summary_${group}.pdf" "matrix_scheme_summary_${group}_page"

    """
}

// --- MULTIQC_REPORT
process MULTIQC_REPORT {
    tag "${group}"

    input:
    tuple val(group),
          path(funnel_png, stageAs: 'variant_filter_funnel.png'),
          path(matrix_pages),
          path(redundancy_png, stageAs: 'hqstat_qc_redundancy.png'),
          path(vaf_outputs),
          path(sparsity_pages),
          path(qc_influence_png, stageAs: 'qc_filter_influence.png'),
          path(qc_threshold_pngs),
          path(first_round_hexbin_png)
    path(ado_png, stageAs: 'ado_germline_comparison.png')

    output:
    tuple val(group), path("multiqc_report_${group}.html"),       emit: report
    tuple val(group), path("multiqc_report_${group}_data"),       emit: data, optional: true

    script:
    """
    set -euo pipefail
    mkdir -p mqc_in

    # Cohort metrics table is intentionally NOT rendered in the report (the
    # cohort_metrics_${group}.tsv/.json are still published by COHORT_METRICS_SUMMARY).

    # ── images (only staged if present) ──────────────────────────────────────────
    [ -s "variant_filter_funnel.png" ] && cp variant_filter_funnel.png mqc_in/01_variant_filter_funnel_mqc.png || true

    i=1; for p in \$(ls -v matrix_scheme_summary_${group}_page-*.png 2>/dev/null); do
        cp "\$p" "mqc_in/02\${i}_matrix_scheme_page\${i}_mqc.png"; i=\$((i + 1)); done

    [ -s "ado_germline_comparison.png" ] && cp ado_germline_comparison.png mqc_in/03_ado_germline_comparison_mqc.png || true

    [ -s "${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png" ]      && cp "${group}_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png"      mqc_in/04_hexbin_first_round_mqc.png  || true
    [ -s "${group}_hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.png" ]  && cp "${group}_hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.png"  mqc_in/05_hexbin_second_round_mqc.png || true
    [ -s "${group}_threshold_exploration_filtered_shared.png" ]               && cp "${group}_threshold_exploration_filtered_shared.png"               mqc_in/06_threshold_shared_mqc.png    || true
    [ -s "${group}_threshold_exploration_filtered_singleton.png" ]            && cp "${group}_threshold_exploration_filtered_singleton.png"            mqc_in/07_threshold_singleton_mqc.png || true

    # redundancy heatmap moved to AFTER the threshold-exploration sections (was 03)
    [ -s "hqstat_qc_redundancy.png" ] && cp hqstat_qc_redundancy.png mqc_in/08_hqstat_qc_redundancy_mqc.png || true

    i=1; for p in \$(ls -v pileup_sparsity_hexbin_${group}_page-*.png 2>/dev/null); do
        cp "\$p" "mqc_in/09\${i}_sparsity_page\${i}_mqc.png"; i=\$((i + 1)); done

    # ── 10-13: artifact-QC filter calibration (AS / PropClipped / BPPos) ──────────
    [ -s "qc_filter_influence.png" ]                        && cp qc_filter_influence.png                        mqc_in/10_qc_filter_influence_mqc.png        || true
    [ -s "qc_threshold_as_propclipped_${group}.png" ]       && cp "qc_threshold_as_propclipped_${group}.png"     mqc_in/11_qc_threshold_as_propclipped_mqc.png || true
    [ -s "qc_threshold_bppos_${group}.png" ]                && cp "qc_threshold_bppos_${group}.png"              mqc_in/12_qc_threshold_bppos_mqc.png         || true
    [ -s "qc_threshold_passrates_${group}.png" ]            && cp "qc_threshold_passrates_${group}.png"          mqc_in/13_qc_threshold_passrates_mqc.png     || true

    echo "[multiqc_report] staged \$(ls mqc_in | wc -l) custom-content files"

    # ── config + run ─────────────────────────────────────────────────────────────
    cat > multiqc_config.yaml <<CFG
title: "Somatic SNP/INDEL Filter Report"
subtitle: "Group: ${group}"
intro_text: "Per-group variant-filtering QC: the funnel, matrix-scheme summary, ADO comparison, Rho-vs-qval hexbins (first + second Sequoia rounds), threshold exploration, HQStat-vs-QC redundancy, pileup VAF sparsity, and artifact-QC filter calibration (influence + AS/PropClipped/BPPos threshold distributions)."
custom_logo_title: "${group}"
show_analysis_paths: false
CFG

    multiqc mqc_in -c multiqc_config.yaml -n multiqc_report_${group}.html -f

    echo "[multiqc_report] Done: multiqc_report_${group}.html"
    """
}


// ============================================================================
// JUNE 2026 SOMATIC UPDATE — MODIFIED MODULES (replaced bodies/signatures)
// ============================================================================

// --- PREPROCESS_VCF
process PREPROCESS_VCF {
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(input_vcf), val(group)
    path(reference)
    val(model_vcf)

  
    output:
    tuple val(sample_name), path("${sample_name}_noref_norm.vcf.gz*"), val(group), emit: vcf
    tuple val(sample_name), path("df_query_${sample_name}.tsv"), val(group), emit: query_table
    tuple val(sample_name), path("df_gt_${sample_name}.tsv"), val(group), emit: df_gt

    
    script:
    """

    if [ "${model_vcf}" = "deepvariant" ]; then

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa | bcftools norm --threads ${task.cpus} -d exact | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz
    else

        bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' ${input_vcf[0]} | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa  | bcftools view --threads ${task.cpus} -i 'GT[*]="alt"' -Oz -o temp.vcf.gz

    fi

    bcftools index --threads ${task.cpus} -t temp.vcf.gz

    echo -e "${sample_name}" > noms.txt;

    bcftools reheader --threads ${task.cpus} -s noms.txt temp.vcf.gz | bcftools view --threads ${task.cpus} -Oz -o ${sample_name}_noref_norm.vcf.gz

    bcftools index --threads ${task.cpus} -t ${sample_name}_noref_norm.vcf.gz

    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%DP\\t%GT]\\n' ${sample_name}_noref_norm.vcf.gz | tail -n+2 >> df_nv_nr.tsv

    cat df_nv_nr.tsv  | awk -v OFS="\\t" '{gsub(".*,","",\$5);print \$1"_"\$2"_"\$3"_"\$4,\$5,\$6}'  > df_query_${sample_name}.tsv

    cat df_nv_nr.tsv  | awk -v OFS="\\t" '{print \$1"_"\$2"_"\$3"_"\$4,\$7}' > df_gt_${sample_name}.tsv
    """
}

// --- MERGE_PROCESSED_VCF
process MERGE_PROCESSED_VCF {
    tag "${group}"
    
    input:
    tuple val(group), path(input_vcfs)
    path(reference)

  
    output:
    tuple val(group), path("merged_group_${group}.vcf.gz*"),emit: merged_vcf
    tuple val(group), path("df_nv_group_${group}.tsv"), emit: df_nv


    script:
    """

    echo -e "Listing files to merge ...";

    find -L . -type f -name "*.vcf.gz" | sort > list_vcf_raw.txt

    cat list_vcf_raw.txt;

    echo -e "Embedding the per-cell pass/reject verdict (FILTER) as a per-sample FORMAT/MLV tag ...";

    # MLV is the per-cell variant-caller verdict carried through the merge as a FORMAT field
    # so that downstream we can count, per variant, in how many cells the variant PASSed vs
    # was rejected. 0 = accepted, 1 = rejected, in that cell. The verdict is caller-agnostic:
    #   - a PASS verdict is FILTER == "PASS" or FILTER == "." (DNAscope writes "." for passing;
    #     DeepVariant writes "PASS")  -> MLV = 0
    #   - any other FILTER value (DNAscope "MLrejected"; DeepVariant "RefCall"/"LowQual"/"NoCall")
    #     is a rejection -> MLV = 1
    # Downstream a variant is kept when it PASSed in >= ml_min_pass_cells cells, i.e. with the
    # default of 1, when at least one cell PASSed it.
    printf '##FORMAT=<ID=MLV,Number=1,Type=Integer,Description="Per-cell variant-caller verdict: 0=accepted (FILTER PASS or .), 1=rejected (any other FILTER, e.g. MLrejected/RefCall/LowQual/NoCall)">\\n' > mlv.hdr

    mkdir -p tagged
    : > list_vcf.txt
    while read -r f; do
        bn=\$(basename "\$f" .vcf.gz)

        # Build a per-cell annotation table mapping each (CHROM,POS,REF,ALT) to its verdict.
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\n' "\$f" \\
            | awk -v OFS='\\t' '{print \$1,\$2,\$3,\$4,((\$5=="PASS" || \$5==".")?0:1)}' \\
            | bgzip > "tagged/\${bn}.mlv.tsv.gz"
        tabix -s1 -b2 -e2 "tagged/\${bn}.mlv.tsv.gz"

        # Annotate adds FORMAT/MLV to the single sample in this per-cell VCF.
        bcftools annotate --threads ${task.cpus} \\
            -a "tagged/\${bn}.mlv.tsv.gz" -h mlv.hdr -c CHROM,POS,REF,ALT,FMT/MLV \\
            "\$f" -Oz -o "tagged/\${bn}.tagged.vcf.gz"
        bcftools index --threads ${task.cpus} -t "tagged/\${bn}.tagged.vcf.gz"

        echo "tagged/\${bn}.tagged.vcf.gz" >> list_vcf.txt
    done < list_vcf_raw.txt

    cat list_vcf.txt;

    echo -e "Merging VCF files ...";

    bcftools merge --threads ${task.cpus} -m none --file-list list_vcf.txt -g ${reference}/genome.fa | bcftools norm --threads ${task.cpus} -m -any --check-ref s -f ${reference}/genome.fa -Oz -o merged_group_${group}.vcf.gz

    echo -e "Indexing merged VCF ...";

    bcftools index --threads ${task.cpus} -t merged_group_${group}.vcf.gz
    
    echo -e "Creating table for NV ...";
    
    bcftools query -l merged_group_${group}.vcf.gz  > columns.txt

    cat columns.txt | tr "\\n" "\\t" | sed "s|\\t\$|\\n|" | sed "s|^|\\t|" > df_nv_group_${group}.tsv

    bcftools query --print-header -f '%CHROM\\_%POS\\_%REF\\_%ALT[\\t%AD]\\n' merged_group_${group}.vcf.gz | tail -n+2 >> df_nv_group_${group}.tsv

    """
}

// --- CUSTOM_BAM_GROUP_PILEUP
process CUSTOM_BAM_GROUP_PILEUP {
    tag "${sample_name}_${chr}"
    
    input:
    tuple val(group), val(sample_name), path(input_bam), path(list_pos), val(chr)
    path(reference)

  
    output:  
    tuple val(group), val(chr), val(sample_name), path("pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt"), emit:pileup
    tuple val(group), val(chr), val(sample_name), path("${sample_name}_df_nr_group${group}_chr_${chr}.tsv"), emit:df_nr

    
    script:
    
    """

    echo -e "Subsetting regions of bam ...";

    date;

    cat ${list_pos} | grep "^${chr}[[:space:]]" | awk -v OFS="\\t" '{print \$1,\$2-1,\$2}' > bed.txt

    # NOTE: `samtools view --regions-file` emits records in BED-row order, not
    # coordinate order, so the resulting BAM is not coordinate-sorted even when
    # the source CRAM/BAM is. That breaks `samtools index` with errors like
    # "Unsorted positions on sequence #N: <pos1> followed by <pos2>". Pipe the
    # subset through `samtools sort` before indexing.
    samtools view --threads $task.cpus --regions-file bed.txt ${input_bam[0]} -O bam -u -o - \\
        | samtools sort --threads $task.cpus -O bam -o regions_subset.bam -

    samtools index --threads $task.cpus regions_subset.bam

    echo -e "Extracting flags/cigar form subsetted bam ...";

    date;

    samtools view --threads ${task.cpus} regions_subset.bam | cut -f1,2,6 | awk -v OFS="\\t" '{print \$1"_"\$2,\$3}' > df_flags_cigars.tsv

    echo -e "Constructing pileup ...";

    date;

    samtools mpileup regions_subset.bam -f ${reference}/genome.fa --output-BP-5 --output-QNAME --disable-overlap-removal --count-orphans --min-BQ 0 --min-MQ 0 --output-MQ --output-extra FLAG,AS -o pileup_temp.txt -l ${list_pos}
    
    echo -e "Adding cigar information to pileup ...";

    date;

    awk -v OFS="\\t" 'NR == FNR {  a[\$1]=\$2; next }{n=split(\$8,ids,",");split(\$9,flags,",");outstring=a[ids[1]"_"flags[1]];for(i=2; i<=n; i++){cigar_string=a[ids[i]"_"flags[i]];outstring=outstring","cigar_string;}print \$0,outstring;}' df_flags_cigars.tsv pileup_temp.txt > pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt

    echo -e "Creating df NR ...";

    date;
    cat pileup_mq0_bq0_group${group}_${sample_name}_${chr}.txt | cut -f1,2,4 | awk -v OFS="\\t" '{print \$1"_"\$2,\$3}' > ${sample_name}_df_nr_group${group}_chr_${chr}.tsv

    """
}

// --- CREATE_TAB_NVNR
process CREATE_TAB_NVNR {
    tag "${group}_${chr}"
    
    input:
    tuple val(group),val(chr), path(df_nv), path(df_nr_files)
  
    output:
    tuple val(group), val(chr), path("mat_nv_group_${group}_chr_${chr}.tsv"), path("mat_nr_group_${group}_chr_${chr}.tsv")

    script:
    """

    echo -e "Listing content to process ...";

    ls;

    echo -e "Cleaning NV ...";

    head -n1 ${df_nv} > header_nv.tsv

    cat ${df_nv} | grep "^${chr}_" | tail -n+2 | awk -v OFS="\\t" '{outstring=\$1;for(i=2; i<=NF; i=i+1){gsub (/.*,/,"",\$i);gsub(/\\./,"0",\$i);outstring=outstring"\\t"\$i;}print outstring;}' | sort -u -k1,1  > df_nv_clean.tsv

    cat header_nv.tsv df_nv_clean.tsv >  mat_nv_group_${group}_chr_${chr}.tsv

    echo -e "Looping of NR files and placing in same order as NV ...";

    cat mat_nv_group_${group}_chr_${chr}.tsv | tail -n+2 | cut -f1 > order_variants_ids_nv.txt

    cat order_variants_ids_nv.txt | cut -d "_" -f1,2 > order_pos_nv.txt

    cat mat_nv_group_${group}_chr_${chr}.tsv| head -n1 | tr "\\t" "\\n" | tail -n+2 > order_samples_nv.txt

    cp order_variants_ids_nv.txt df_nr_clean.tsv

    cat order_samples_nv.txt | while read nom;
    do

        touch temp_nr.tsv;

        ls \${nom}_df_nr_group*.tsv  | sort -V  | while read nr; do cat \${nr} >> temp_nr.tsv;done

        echo -e "Pocessing \${nom}";

        if [ -s temp_nr.tsv ]; then
            awk -v OFS="\\t" 'NR == FNR { a[\$1]=\$2; next }{if( \$0 in a ){print a[\$0]}else{print 0}}' temp_nr.tsv order_pos_nv.txt > temp.txt
        else
            awk '{print 0}' order_pos_nv.txt > temp.txt
        fi

        paste -d "\\t" df_nr_clean.tsv temp.txt > temp_paste.txt

        mv temp_paste.txt df_nr_clean.tsv

        rm temp.txt

        rm temp_nr.tsv;
    
    done

    cat header_nv.tsv df_nr_clean.tsv >  mat_nr_group_${group}_chr_${chr}.tsv

    wc -l mat*.tsv

    nv_rows=\$(tail -n+2 mat_nv_group_${group}_chr_${chr}.tsv | wc -l)
    nr_rows=\$(tail -n+2 mat_nr_group_${group}_chr_${chr}.tsv | wc -l)
    echo "NV data rows: \${nv_rows} | NR data rows: \${nr_rows}"
    if [ "\${nv_rows}" -ne "\${nr_rows}" ]; then
        echo "ERROR: NR row count (\${nr_rows}) does not match NV row count (\${nv_rows})" >&2
        exit 1
    fi

    """
}

// --- CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR
process CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR {
    tag "${group}"
    
    input:
    tuple val(group), path(files), path(priority_variants)
    val(cutoff_binomial)
    val(cutoff_beta)
  
    output:
    tuple val(group), path("res_${group}_binomial_betabinomial.tsv"), emit: res_df
    tuple val(group), path("chosen_variants_${group}.txt"), emit: chosen_variants

    script:
    """

    echo -e "Concatenating tables ...";
    cat filtered*.txt | grep "Depth_filter" | sed "s|^|VariantId\\t|" | sort -u > res.tsv

    ls filtered*.txt | sort -V | while read file; do cat \${file} | tail -n+2 >> res.tsv  ;done

    echo -e "Launching filtering with binomial cutoff: ${cutoff_binomial} and beta-binomial cutoff: ${cutoff_beta}";

    Rscript /usr/local/bin/rscript_0.filter_binom_betabinom_tables.R res.tsv ${cutoff_binomial} ${cutoff_beta}

    mv df_verdict.txt res_${group}_binomial_betabinomial.tsv

    mv chosen_variants.txt original_chosen_${group}.txt

    # Extract rescue IDs from priority_variants TSV (col 1 = VariantId, skip header)
    if [ -s "${priority_variants}" ]; then
      awk 'NR>1{print \$1}' ${priority_variants} > rescue_ids.txt
    else
      touch rescue_ids.txt
    fi

    # Union: original statistical passes + rescued priority variants
    cat original_chosen_${group}.txt rescue_ids.txt | sort -u > chosen_variants_${group}.txt

    echo -e "Statistical filter retained \$(wc -l < original_chosen_${group}.txt) variants"
    echo -e "Rescued \$(comm -13 <(sort original_chosen_${group}.txt) <(sort rescue_ids.txt) | wc -l) additional priority variants"

    # Annotate provenance: Pass (original), Rescued (added by priority), Fail (removed)
    awk -v OFS="\\t" \\
        'FILENAME==ARGV[1] { orig[\$0]=1; next }
         FILENAME==ARGV[2] { rescue[\$0]=1; next }
         FNR==1 { print \$0, "BinomialBetabinomialFilter"; next }
         {
           if      (\$1 in orig)   status="Pass";
           else if (\$1 in rescue) status="Rescued";
           else                   status="Fail";
           print \$0, status
         }' original_chosen_${group}.txt rescue_ids.txt res_${group}_binomial_betabinomial.tsv > res_annotated.tsv

    mv res_annotated.tsv res_${group}_binomial_betabinomial.tsv

    """
    
}

// --- CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR
process CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR {
    tag "${sample_name}_${chr}"
    

    input:
    tuple val(group), val(chr), val(sample_name), path(pileup_cigars_file), path(chosen_variants)
    val(threshold_mq)
    val(threshold_bq)
    val(threshold_bp)
    val(num_lines_read_pileup)
    val(read_length)
  
    output:
    tuple val(group),val(chr), path("res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv")

    script:
    """

    echo -e "Subsetting variants ...";

    cat ${chosen_variants} | grep "^${chr}_" > cvariants.txt || true

    cat cvariants.txt | cut -d "_" -f2 | sort -uV > cpositions.txt
    
    awk -v FS="\\t" -v OFS="\\t" 'NR == FNR {  a[\$1]=1; next }{if( \$2 in a ){print \$0}}' cpositions.txt ${pileup_cigars_file} > pileup_subset.tsv

    echo -e "Launching QC script ...";

    if [ -s "pileup_subset.tsv" ]; then
        Rscript /usr/local/bin/rscript_1.samplelevel_process_pileup_sample_cigars.R pileup_subset.tsv ${sample_name} ${threshold_mq} ${threshold_bq} ${threshold_bp} ${num_lines_read_pileup} ${read_length} || true
        if [ -f "res.tsv" ]; then
            mv res.tsv res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv
        else
            echo -e "R script produced no output (likely all READ_BASES empty), writing empty result ...";
            touch res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv
        fi
    else
        echo -e "Pileup subset file is empty, skipping QC script ...";
        touch res_grouplevel_pileup_group_${group}_sample_${sample_name}_chr_${chr}.tsv;
    fi
    
    """
    
}

// --- CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES
process CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES {
    tag "${group}_${chr}"
    
    input:
    tuple val(group),val(chr), path(res_tables), path(chosen_variants), path(priority_variants)
    val(threshold_as)
    val(threshold_clipped)
    val(threshold_prop_bp_under)  
    val(threshold_prop_bp_upper)  
    val(threshold_sd_indiv)
    val(threshold_mad_indiv)
    val(threshold_sd_both)
    val(threshold_mad_both)
    val(threshold_sd_extreme)
    val(threshold_mad_extreme)
    val(disable_qc)
    val(disable_bppos)


    output:
    tuple val(group), val(chr), path("Mat_NV_${group}_${chr}.tsv"), path("Mat_NR_${group}_${chr}.tsv"), emit: tabs
    tuple val(group), val(chr), path("Mat_NV_all_${group}_${chr}.tsv"), path("Mat_NR_all_${group}_${chr}.tsv"), emit: tabs_all
    tuple val(group), val(chr), path("df_passed*"), emit:df_pass
    tuple val(group), val(chr), path("res_pileup_all_group_${group}_${chr}.tsv"), emit:pileup

    script:
    
    """



    echo -e "Concatenating tables ...";

    find . -name "res_grouplevel_pileup*" | sort -V > list_files.txt

    echo -e "SampleId\\tVariantId\\tCHROM\\tPOS\\tREF\\tALT\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_F\\tNUM_FRAGMENTS_ALLQ_MQ_BQ_R\\tNUM_FRAGMENTS_HQ_MQ_BQ_F\\tNUM_FRAGMENTS_HQ_MQ_BQ_R\\tMEDIAN_AS_VARIANT_READS\\tPROP_BASES_CLIPPED\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F\\tPROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F\\tSD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R\\tMAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F\\tMAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R\\tNUM_FRAGMENTS_ALLQ_POSITION\\tNUM_FRAGMENTS_HQ_POSITION\\tPROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F\\tPROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R" > res_end.tsv;

    cat list_files.txt | while read mfile; do cat \${mfile} >> res_end.tsv;done

    head -n1 res_end.tsv > df_raw_variants.tsv
    
    cat res_end.tsv | grep -v VariantId | grep -Pv "\\tREF\\t" >> df_raw_variants.tsv || true

    echo -e "Filtering ... ";

    Rscript /usr/local/bin/rscript_2.create_tabnr_tabnv.R ${threshold_as} ${threshold_clipped} ${threshold_prop_bp_under} ${threshold_prop_bp_upper} ${threshold_sd_indiv} ${threshold_mad_indiv} ${threshold_sd_both} ${threshold_mad_both} ${threshold_sd_extreme} ${threshold_mad_extreme} ${disable_qc} ${disable_bppos}

    echo -e "Subsetting matrices ... ";

    head -n1 Tab_NR.tsv > Mat_NR_${group}_${chr}.tsv

    tail -n +2 Tab_NR.tsv > body.tsv

    cat ${chosen_variants} | grep "^${chr}_" > mvariants.txt || true

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' mvariants.txt body.tsv >> Mat_NR_${group}_${chr}.tsv

    # Rescue priority variants into Mat_NR: pull from Tab_NR_all rows in mvariants that QC filtering removed
    if [ -s "${priority_variants}" ]; then
      awk 'NR>1{print \$1}' ${priority_variants} | grep "^${chr}_" > rescue_chr.txt || true
      if [ -s rescue_chr.txt ]; then
        tail -n+2 Mat_NR_${group}_${chr}.tsv | cut -f1 | sort > already_in_nr.txt
        comm -23 <(sort rescue_chr.txt) <(sort already_in_nr.txt) > to_rescue_nr.txt
        awk 'NR==FNR{r[\$0]=1;next} \$0 in r' to_rescue_nr.txt mvariants.txt > rescue_nr_filtered.txt
        if [ -s rescue_nr_filtered.txt ]; then
          awk -v OFS="\\t" -v FS="\\t" 'NR==FNR{a[\$0]=1;next} FNR>1 && (\$1 in a){print \$0}' rescue_nr_filtered.txt Tab_NR_all.tsv >> Mat_NR_${group}_${chr}.tsv
          echo "Rescued \$(wc -l < rescue_nr_filtered.txt) priority variants into Mat_NR_${group}_${chr}"
        fi
      fi
    fi

    cat Mat_NR_${group}_${chr}.tsv | tail -n +2 | cut -f1 > vfound.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(!(\$1 in a)){ print\$0}}' vfound.txt mvariants.txt > df_passed_PRESENTVCF_NOTINBAM_${group}_${chr}.tsv

    cat body.tsv  | cut -f1 > foundbam.txt

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(!(\$1 in a)){ print \$0}}' mvariants.txt foundbam.txt > df_passed_PRESENTBAM_NOTINVCF_${group}_${chr}.tsv

    head -n1 Tab_NV.tsv > Mat_NV_${group}_${chr}.tsv

    tail -n +2 Tab_NV.tsv > body.tsv

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' ${chosen_variants} body.tsv >> Mat_NV_${group}_${chr}.tsv

    # Rescue priority variants into Mat_NV: pull from Tab_NV_all rows in chosen_variants that QC filtering removed
    if [ -s "${priority_variants}" ]; then
      awk 'NR>1{print \$1}' ${priority_variants} > rescue_ids.txt
      if [ -s rescue_ids.txt ]; then
        tail -n+2 Mat_NV_${group}_${chr}.tsv | cut -f1 | sort > already_in_nv.txt
        comm -23 <(sort rescue_ids.txt) <(sort already_in_nv.txt) > to_rescue_nv.txt
        awk 'NR==FNR{r[\$0]=1;next} \$0 in r' to_rescue_nv.txt ${chosen_variants} > rescue_nv_filtered.txt
        if [ -s rescue_nv_filtered.txt ]; then
          awk -v OFS="\\t" -v FS="\\t" 'NR==FNR{a[\$0]=1;next} FNR>1 && (\$1 in a){print \$0}' rescue_nv_filtered.txt Tab_NV_all.tsv >> Mat_NV_${group}_${chr}.tsv
          echo "Rescued \$(wc -l < rescue_nv_filtered.txt) priority variants into Mat_NV_${group}_${chr}"
        fi
      fi
    fi

    echo -e "Subsetting all-variants matrices ...";

    head -n1 Tab_NR_all.tsv > Mat_NR_all_${group}_${chr}.tsv

    tail -n +2 Tab_NR_all.tsv > body_all.tsv

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' mvariants.txt body_all.tsv >> Mat_NR_all_${group}_${chr}.tsv

    head -n1 Tab_NV_all.tsv > Mat_NV_all_${group}_${chr}.tsv

    tail -n +2 Tab_NV_all.tsv > body_nv_all.tsv

    awk -v OFS="\\t" -v FS="\\t" 'NR == FNR {  a[\$0]; next }{if(\$1 in a){ print \$0}}' ${chosen_variants} body_nv_all.tsv >> Mat_NV_all_${group}_${chr}.tsv

    echo -e "Renaming files to output ...";
    mv df_passed_AS.tsv df_passed_AS_${group}_${chr}.tsv

    mv df_passed_propclipped.tsv df_passed_propclipped_${group}_${chr}.tsv

    mv df_passed_BPPOS.tsv df_passed_BPPOS_${group}_${chr}.tsv

    cat Mat_NV_${group}_${chr}.tsv | tail -n +2 | cut -f1 > df_passed_DEPTH_${group}_${chr}.tsv

    mv res_end.tsv res_pileup_all_group_${group}_${chr}.tsv

    echo -e "Annotating pileup table with per-filter status ...";

    awk -v OFS="\\t" '
      FILENAME ~ /df_passed_AS/                  && NF>=3 { as_pass[\$2 SUBSEP \$3]=1;      next }
      FILENAME ~ /df_passed_propclipped/         && NF>=3 { clip_pass[\$2 SUBSEP \$3]=1;    next }
      FILENAME ~ /df_passed_BPPOS/               && NF>=3 { bppos_pass[\$2 SUBSEP \$3]=1;   next }
      FILENAME ~ /df_passed_DEPTH/               && NF>=1 { depth_pass[\$1]=1;             next }
      FILENAME ~ /df_passed_PRESENTVCF_NOTINBAM/ && NF>=1 { vcf_notbam[\$1]=1;             next }
      FILENAME ~ /df_passed_PRESENTBAM_NOTINVCF/ && NF>=1 { bam_notvcf[\$1]=1;             next }
      FNR==1 {
        print \$0, "AS_Filter", "PropClipped_Filter", "BPPos_Filter",
                   "Depth_Filter", "PresentVCF_NotInBAM", "PresentBAM_NotInVCF", "Verdict"
        next
      }
      {
        key = \$1 SUBSEP \$2; vid = \$2
        print \$0,
          ((key in as_pass)      ? "Pass" : "Fail"),
          ((key in clip_pass)    ? "Pass" : "Fail"),
          ((key in bppos_pass)   ? "Pass" : "Fail"),
          ((vid in depth_pass)   ? "Pass" : "Fail"),
          ((vid in vcf_notbam)   ? "Yes"  : "No"),
          ((vid in bam_notvcf)   ? "Yes"  : "No"),
          ((vid in depth_pass)   ? "Pass" : "Fail")
      }
    ' df_passed_AS_${group}_${chr}.tsv \
      df_passed_propclipped_${group}_${chr}.tsv \
      df_passed_BPPOS_${group}_${chr}.tsv \
      df_passed_DEPTH_${group}_${chr}.tsv \
      df_passed_PRESENTVCF_NOTINBAM_${group}_${chr}.tsv \
      df_passed_PRESENTBAM_NOTINVCF_${group}_${chr}.tsv \
      res_pileup_all_group_${group}_${chr}.tsv > res_pileup_annotated.tsv

    mv res_pileup_annotated.tsv res_pileup_all_group_${group}_${chr}.tsv

    """
}

// --- CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS
process CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS {
    tag "${group}"
    
    input:
    tuple val(group), path(res_tables_nv),  path(res_tables_nr)


    output:
    tuple val(group), path("Mat_NV_${group}.tsv"), path("Mat_NR_${group}.tsv"), emit: tabs

    script:

    """

    for f in Mat_NR*; do [ \$(wc -l < "\$f") -gt 1 ] && echo "\$f"; done | sort -V > list_nr.txt

    for f in Mat_NV*; do [ \$(wc -l < "\$f") -gt 1 ] && echo "\$f"; done | sort -V > list_nv.txt

    ls Mat_N* | while read file; do head -n1 \${file} | tr "\\t" "\\n" | tail -n+2;done | sort -uV > all_samples.txt

    Rscript /usr/local/bin/rscript_3.create_tabnr_tabnv_group.R

    mv Tab_NR.tsv Mat_NR_${group}.tsv

    mv Tab_NV.tsv Mat_NV_${group}.tsv


    """
}

// --- CUSTOM_VARIANT_FILTER_PROVENANCE
// Build a master variant × filter-status table by joining all per-stage provenance outputs.
// Also generates a filter tracking table showing how many variants are retained at each stage.
// Join order: bulk → binom → vep → pileup → sequoia.
// All columns from every source table are kept, prefixed with the source name.
// Pileup is sample-level: numeric metrics are averaged, filter cols become PassCounts,
// SampleId becomes NumSamples, and per-variant-constant fields take their first value.
process CUSTOM_VARIANT_FILTER_PROVENANCE {
    tag "${group}"

    input:
    // prefilter_table = variant_prefilter_table_${group}.tsv (Pileup_* + SecondRunSequoia_*,
    //   consumed by master_table.R instead of recomputing them).
    // vaf_split_outputs = VAF_SPLIT_VARIANTS_HEXBIN.out.all_outputs file set; master_table.R
    //   globs the singleton/shared/retained lists + cascade binary matrices for VAFSplit_* columns.
    // (all_variants, the sequoia df, and the pileup pool are no longer inputs here:
    //  master_table.R reads bulk/binom/vep + prefilter; the sequoia df + pileup pool are
    //  consumed by VARIANT_PREFILTER_TABLE; the pileup pool also by BUILD_FOCAL_PILEUP;
    //  all_variants only by VARIANT_FILTER_FUNNEL.)
    tuple val(group), path(bulk_prov), path(binom_tsv), path(vep_prov), path(tab_nvnr_files), path(df_nv_prebulk), path(priority_variants), path(prefilter_table), path(vaf_split_outputs)

    output:
    tuple val(group), path("variant_master_filter_table_${group}.tsv"),        emit: master_table
    // QC-plot + combined-report emits disabled (scripts commented out in script block):
    //   binom_plots, vaf_hexbin_plots, pileup_metric_plots, pileup_bppos_plots, combined_report
    tuple val(group), path("vcf_annotation_table_${group}.tsv"),               emit: vcf_annotation_table
    // pileup_focal moved to the BUILD_FOCAL_PILEUP module (the streaming focal-pileup awk).
    tuple val(group), path("upstream_filter_per_sample_${group}.tsv"),        emit: upstream_per_sample

    script:
    def g = group
    """
    set -euo pipefail


    Rscript /usr/local/bin/master_table.R "${g}"

    # Join priority_variants annotation columns onto master table.
    # Adds Impact, Consequence, GenesAffected, DataBasesReported, Verdict, PriorityVariant (Yes/No).
    # Variants absent from priority_variants get "." for annotation columns and "No" for PriorityVariant.
    awk -F'\\t' -v OFS='\\t' '
      FILENAME==ARGV[1] {
        if (FNR==1) { ncols=NF; for(i=2;i<=NF;i++) hdr[i]=\$i; next }
        for(i=2;i<=ncols;i++) pv[\$1,i]=\$i
        pv_flag[\$1]="Yes"
        next
      }
      FNR==1 {
        printf "%s", \$0
        for(i=2;i<=ncols;i++) printf "\\t%s", hdr[i]
        print "\\tPriorityVariant"
        next
      }
      {
        printf "%s", \$0
        if (\$1 in pv_flag) {
          for(i=2;i<=ncols;i++) printf "\\t%s", pv[\$1,i]
          print "\\tYes"
        } else {
          for(i=2;i<=ncols;i++) printf "\\t."
          print "\\tNo"
        }
      }
    ' ${priority_variants} variant_master_filter_table_${g}.tsv > master_annotated.tsv
    mv master_annotated.tsv variant_master_filter_table_${g}.tsv

    echo "Master table annotated with \$(awk -F'\\t' 'NR>1 && \$NF==\"Yes\"' variant_master_filter_table_${g}.tsv | wc -l) priority variants"

    # Plotting/post-processing scripts run sequentially. Running them in parallel
    # caused OOM kills (exit 255) because pileup_metric_plots.R and
    # pileup_bppos_plots.R each rbind ~140M-row pileup pools while the master
    # table (16M+ rows) is also being re-read by other scripts in the same R
    # process group. Serialising lets each script start with a clean R heap.
    Rscript /usr/local/bin/upstream_filter_per_sample.R \
        "${g}" \
        variant_master_filter_table_${g}.tsv \
        "${df_nv_prebulk}" \
        mat_nv_group_*.tsv

    # ── QC plots + combined report DISABLED for now ──────────────────────────
    # The four per-variant QC-plot scripts and the pdfunite step are commented out;
    # the module no longer emits binom/vaf_hexbin/pileup_metric/pileup_bppos PDFs or
    # the combined report. Re-enable by uncommenting these and restoring the emits.
    # Rscript /usr/local/bin/binom_ggplots.R variant_master_filter_table_${g}.tsv "${g}" binom_filter_plots_${g}.pdf
    #
    # Rscript /usr/local/bin/vaf_hexbin_plots.R variant_master_filter_table_${g}.tsv "${g}" vaf_hexbin_plots_${g}.pdf
    #
    # Rscript /usr/local/bin/pileup_metric_plots.R "${g}" pileup_metric_plots_${g}.pdf
    #
    # Rscript /usr/local/bin/pileup_bppos_plots.R "${g}" variant_master_filter_table_${g}.tsv pileup_bppos_plots_${g}.pdf
    #
    # pdfunite \
    #     binom_filter_plots_${g}.pdf \
    #     vaf_hexbin_plots_${g}.pdf \
    #     pileup_metric_plots_${g}.pdf \
    #     pileup_bppos_plots_${g}.pdf \
    #     Patient_filter_report_combined_${g}.pdf

    Rscript /usr/local/bin/vcf_annotation_table.R "${g}"

    """
}

// --- SEQUOIA_SECOND_FILTER
process SEQUOIA_SECOND_FILTER {
    tag "${group}_${chr}"

    input:
    tuple val(group), val(chr), path(mat_nv), path(mat_nr)
    path(reference)
    val(cutoff_binomial)
    val(cutoff_rho_snp)
    val(cutoff_rho_indel)
    val(min_cov)
    val(max_cov)
    val(gender)
    val(beta_binom_shared)

    output:
    tuple val(group), val(chr), path("*_filtering_all.txt"), emit: df_filter

    script:
    """

    echo -e "Subsetting matrices to chromosome ${chr} ..."

    awk -v CHR="${chr}" 'NR==1 || \$1 ~ "^"CHR"_"' ${mat_nv} > mat_nv_chr.tsv
    awk -v CHR="${chr}" 'NR==1 || \$1 ~ "^"CHR"_"' ${mat_nr} > mat_nr_chr.tsv

    data_rows=\$(tail -n+2 mat_nv_chr.tsv | wc -l)
    echo "Chromosome ${chr}: \${data_rows} variant rows"


    if [ "\${data_rows}" -gt 0 ]; then

        echo -e "Raw matrices number of lines ...";
        wc -l mat_nv_chr.tsv
        wc -l mat_nr_chr.tsv

        Rscript /usr/local/bin/rscript_4.sequoia_second_pass_filter.R --genomeFile ${reference}/genome.fa -v mat_nv_chr.tsv -r mat_nr_chr.tsv --mpboot_path /usr/local/bin/ -n $task.cpus --snv_rho ${cutoff_rho_snp} --indel_rho ${cutoff_rho_indel} --germline_cutoff ${cutoff_binomial} --min_cov ${min_cov} --max_cov ${max_cov} --gender ${gender} -b ${beta_binom_shared}

        ls Patient* | while read file;
        do

            name=`echo \${file} | sed 's/Patient/Sequoia_group_${group}_chr_${chr}_bino${cutoff_binomial}_rhosnp${cutoff_rho_snp}_rhoindel${cutoff_rho_indel}_mincov${min_cov}_maxcov${max_cov}/'`;

            mv \${file} \${name};

        done

        echo -e "Annotating filtering table with SecondPassFilter column ..."

        filt_file=\$(ls *_filtering_all.txt)
        nr_file=\$(ls *_NR_filtered_all.txt)

        awk 'NR==FNR { if (FNR>1) pass[\$1]=1; next }
             FNR==1  { print \$0, "SecondPassFilter"; next }
             { print \$0, ((\$1 in pass) ? "Pass" : "Fail") }
        ' "\${nr_file}" "\${filt_file}" > filtering_annotated.txt

        mv filtering_annotated.txt "\${filt_file}"

    else

        echo "No variants for chromosome ${chr} — creating empty output"
        echo -e "variantID\tSecondPassFilter" > Sequoia_group_${group}_chr_${chr}_bino${cutoff_binomial}_rhosnp${cutoff_rho_snp}_rhoindel${cutoff_rho_indel}_mincov${min_cov}_maxcov${max_cov}_filtering_all.txt

    fi

    """

}

// --- SEQUOIA_SECOND_FILTER_MERGE
process SEQUOIA_SECOND_FILTER_MERGE {
    tag "${group}"

    input:
    tuple val(group), path(filtering_files)

    output:
    tuple val(group), path("Sequoia_group_${group}_filtering_all.txt"), emit: df_filter

    script:
    """
    # Take header from the first file; all per-chromosome files share the same header
    first_file=\$(ls *_chr_*_filtering_all.txt | sort -V | head -1)
    head -n1 "\${first_file}" > Sequoia_group_${group}_filtering_all.txt

    # Append data rows from all per-chromosome files in sorted order
    for f in \$(ls *_chr_*_filtering_all.txt | sort -V); do
        tail -n+2 "\${f}" >> Sequoia_group_${group}_filtering_all.txt
    done
    """
}

// ---------------------------------------------------------------------------
// EMIT_LINEAGE_INPUTS_MANIFEST
// Emit the basej-lineage handoff manifest(s) for the downstream connector in one
// pass (pure path templating — no file staging):
//   * index/lineage_inputs.csv       — "group,param,path" (one row-set per group;
//                                       multi-cohort safe, no filename collision;
//                                       row-extensible for future inputs)
//   * index/lineage_vcfs_<group>.csv — per-group "biosampleName,vcf,vcf_index"
//                                       (vcf_somatic is a shared flat folder, so a
//                                       per-group input_csv scopes each lineage run
//                                       to its own cells; referenced by the manifest)
// The connector reads these paths verbatim — no base-prefix assembly, no guessing
// which filtered matrix is canonical (we encode the choice here). NR/NV/binary point
// at the COLLECT_DOWNSTREAM_ARTIFACTS bundle's PlusMandatoryNonEmpty scheme (see the
// script body for the rationale — ForPhylogeny starves the lineage heatmap step).
//   group_samples = list of [group, [samples]] (only groups that produced matrices)
//   out_base      = <outputDir>/workflow_outputs/<workspace>/<workflow_id>
// Container + resources are configured in nextflow.config under withName.
// ---------------------------------------------------------------------------
process EMIT_LINEAGE_INPUTS_MANIFEST {
    tag "manifest"

    input:
    val(group_samples)
    val(out_base)

    output:
    path("lineage_inputs.csv"),  emit: manifest
    path("lineage_vcfs_*.csv"),  emit: vcf_csvs, optional: true

    script:
    def sa      = "${out_base}/secondary_analyses"
    // NR/NV/binary matrices for lineage are the COLLECT_DOWNSTREAM_ARTIFACTS bundle copies
    // (secondary_analyses/downstream_inputs/downstream_inputs_<group>/), NOT phylogeny_matrices/.
    // We hand off the PlusMandatoryNonEmpty scheme (mandatory variants included, and WITHOUT the
    // phylogeny down-selection). Rationale: the pre-selected ForPhylogeny matrix is placed in full
    // by SEQUOIA_PHYLOGENY, which leaves the downstream variant-placement step with nothing to place
    // and therefore SKIPS the VAF/digital heatmaps. The fuller PlusMandatoryNonEmpty matrix leaves
    // variants for placement, so the heatmaps render (matches Isai's lineage runs).
    def suf     = "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredPlusMandatoryNonEmpty"
    def entries = group_samples instanceof List ? group_samples : [group_samples]
    def mlines  = ["group,param,path"]
    def vcf_blocks = []
    entries.each { e ->
        def g       = e[0]
        def samples = e[1] instanceof List ? e[1] : [e[1]]
        def dsdir   = "${sa}/downstream_inputs/downstream_inputs_${g}"
        mlines << "${g},nr_matrix,${dsdir}/${g}_NR_${suf}.tsv"
        mlines << "${g},nv_matrix,${dsdir}/${g}_NV_${suf}.tsv"
        mlines << "${g},binary_matrix,${dsdir}/${g}_binary_matrix_${suf}.tsv"
        mlines << "${g},mandatory_variants_qc_status,${sa}/mandatory_qc/mandatory_variants_qc_status_${g}.tsv"
        mlines << "${g},input_csv,${out_base}/index/lineage_vcfs_${g}.csv"
        def vlines = ["biosampleName,vcf,vcf_index"]
        samples.each { s ->
            vlines << "${s},${sa}/vcf_somatic/${s}_somatic_annotated.vcf.gz,${sa}/vcf_somatic/${s}_somatic_annotated.vcf.gz.tbi"
        }
        vcf_blocks << "cat > lineage_vcfs_${g}.csv <<'VEOF'\n${vlines.join('\n')}\nVEOF"
    }
    def manifest_body = mlines.join('\n')
    def vcf_cmds      = vcf_blocks.join('\n')
    """
    cat > lineage_inputs.csv <<'MEOF'
${manifest_body}
MEOF
${vcf_cmds}
    """
}


