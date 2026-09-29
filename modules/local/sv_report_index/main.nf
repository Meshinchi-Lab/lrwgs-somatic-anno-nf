process SV_REPORT_INDEX {
    tag "cohort"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/rocker/verse:4.6.1' :
        'docker.io/rocker/verse:4.6.1' }"

    input:
    path annotsv_tsvs
    path knot_htmls,            stageAs: 'knot/*'
    path knot_xls,              stageAs: 'knot/*'   // .xlsm files land in the same dir as the HTML reports
    path vep_html_reports,      stageAs: 'vep/*'
    path igv_html_reports,      stageAs: 'igv/*'
    path vembrane_pass2_tsvs
    path vembrane_pass2_cna_tsvs, stageAs: 'cna_pass2/*'   // CNA pass2 TSVs (FILTER_CNA.out.tsv_pass2); stagedAs avoids glob collision with SV `.filtered.pass2.tsv`
    path vembrane_snv_tsvs,     stageAs: 'snv/*'   // optional; from FILTER_SNV
    path sample_metadata
    path candidate_genes        // optional one-gene-per-line file
    path ensembl_db,            stageAs: 'ensembl_db/*'      // optional EnsDb sqlite (or similar); consumed by the qmd via params$ensembl_db
    path annotsv_annotations,   stageAs: 'annotsv_annotations'  // AnnotSV annotations dir (symlinked); qmd globs Users/**/*.header.tsv to derive the ST regex
    // Per-sample VCF header dictionaries emitted by CREATE_DATA_DICT_{SV,CNA,SNV}.
    // Staged into dedicated subdirs so the qmd can glob each variant type's dicts
    // and union them (identical field_names across samples are collapsed to one
    // definition row). optional:true keeps stub runs and single-branch runs green.
    path data_dict_sv,          stageAs: 'data_dict_sv/*'
    path data_dict_cna,         stageAs: 'data_dict_cna/*'
    path data_dict_snv,         stageAs: 'data_dict_snv/*'
    // Cohort-merge results (TWO-PASS RENDER -- see the workflow call site).
    // `driver_tier` is computed in the Quarto layer, so COHORT_MERGE_SV cannot
    // run until this process has emitted cohort_sv_pass2.tsv once. A single
    // render therefore cannot display its own cohort results. Pass 1 receives
    // `[]` for both and renders a "pending" card; pass 2 receives the merge
    // outputs and supersedes pass 1 in the publish dir.
    path cohort_recurrent_tsv,  stageAs: 'cohort_merge/*'
    path cohort_merged_vcf,     stageAs: 'cohort_merge/*'
    // DISCO on/off for THIS render, set at the call site rather than read from
    // params, because the two passes need different values and a withName
    // selector cannot express it. published report, rendered with DISCO off. `null` falls back to params.disco_enable.
    val disco_enable_in

    output:
    path "index.html",                                           emit: report
    path "sv_report_index.qmd",                                  emit: qmd
    path "*.html",                                               emit: html_plots
    path "*.html.log",                                           emit: log
    path "cohort_sv_pass2.tsv",  optional: true,   emit: tsv_sv
    path "cohort_cna_pass2.tsv", optional: true,   emit: tsv_cna
    path "cohort_snv.tsv", optional: true,   emit: tsv_snv
    // Per-sample ProteinPaint DISCO HTMLs + per-gene Tier-1 SNV lollipop
    // HTMLs. The qmd writes both families under `protein_paint_out/` so a
    // single output declaration covers both. `optional: true` keeps the
    // process green when disco_enable=false or no Tier-1 records survive.
    path "protein_paint_out/*.html", optional: true, emit: protein_paint_html
    // ProteinPaint gene-fusion input (Tier 1 + Tier 2 SVs, svfusion format).
    // Separate declaration because the HTML glob above does not match *.txt.
    path "protein_paint_out/*.svfusion.txt", optional: true, emit: sv_fusion_pp
    // Re-publish the per-variant-type VCF header data dictionaries that
    // were staged into this task as inputs. Making them outputs of
    // SV_REPORT_INDEX causes publishDir to copy them into
    // `results/sv_report/data_dict_{sv,cna,snv}/` alongside index.html,
    // which is the relative path the qmd's "Download full data dictionary"
    // hyperlinks resolve against. optional:true tolerates a disabled
    // CNA / SNV branch (empty staging dir → no glob match).
    path "data_dict_sv/*.data_dict.tsv",  optional: true, emit: data_dict_sv_out
    path "data_dict_cna/*.data_dict.tsv", optional: true, emit: data_dict_cna_out
    path "data_dict_snv/*.data_dict.tsv", optional: true, emit: data_dict_snv_out
    tuple val("${task.process}"), val('quarto'),
          eval("quarto --version"),
          emit: versions_quarto, topic: versions
    tuple val("${task.process}"), val('r-base'),
          eval("Rscript -e 'cat(R.version[[\"major\"]], \".\", R.version[[\"minor\"]], sep=\"\")'"),
          emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Path values are passed to the qmd via `quarto render -P key:value` and
    // consumed as `params$meta_data` / `params$candidate_genes` / `params$ensembl_db`.
    // Empty string when the optional input isn't provided (the qmd checks nzchar() + file.exists() before using them).
    def metadata_arg  = sample_metadata  ? sample_metadata.toString()  : ""
    def candgenes_arg = candidate_genes  ? candidate_genes.toString()  : ""
    def ensdb_arg     = ensembl_db       ? ensembl_db.toString()       : ""
    // ProteinPaint DISCO plot params — nextflow.config defaults are
    // pipeline-wide; profile overrides (test.config, OPBG) can turn the
    // plot off or point to a different host without touching the qmd.
    // Explicit input (see the `disco_enable_in` comment above); null-safe so a
    // caller that does not care still honours params.disco_enable.
    def disco_flag         = disco_enable_in != null ? disco_enable_in : params.disco_enable
    def disco_enable_arg   = disco_flag ? 'true' : 'false'
    // Empty string on pass 1 -> the qmd renders a "pending" note, not the tables.
    def cohort_tsv_arg     = cohort_recurrent_tsv ? cohort_recurrent_tsv.toString() : ""
    def cohort_vcf_arg     = cohort_merged_vcf    ? cohort_merged_vcf.toString()    : ""
    def pp_host_arg        = params.proteinpaint_host   ?: 'http://localhost:3456'
    def pp_genome_arg      = params.proteinpaint_genome ?: 'hg38'
    """
    export HOME="\${PWD}"
    cp ${projectDir}/bin/sv_report_index.qmd .

    quarto render sv_report_index.qmd \\
        --output index.html \\
        --no-cache \\
        -P meta_data:"${metadata_arg}" \\
        -P candidate_genes:"${candgenes_arg}" \\
        -P ensembl_db:"${ensdb_arg}" \\
        -P disco_enable:${disco_enable_arg} \\
        -P proteinpaint_host:"${pp_host_arg}" \\
        -P proteinpaint_genome:"${pp_genome_arg}" \\
        -P cohort_recurrent_tsv:"${cohort_tsv_arg}" \\
        -P cohort_merged_vcf:"${cohort_vcf_arg}" \\
        --log knotannotsv_vep.html.log \\
        --log-level debug
    """

    stub:
    """
    touch index.html
    touch sv_report_index.qmd
    touch knotannotsv_vep.html.log

    # Emit a structurally valid cohort_sv_pass2.tsv so stub runs actually exercise
    # downstream consumers (COHORT_MERGE_SV). The emit is `optional: true`, so without
    # this the channel is empty and the whole cohort merge silently disappears from the
    # DAG -- indistinguishable from it being disabled.
    printf 'SAMPLE\tCHROM\tPOS\tCytoBand\tID\tREF\tALT\tSVTYPE\tDETAILED_TYPE\tSVLEN\tEND\tdriver_tier\n' > cohort_sv_pass2.tsv
    printf 'stub\t21\t33932155\t21q22.12\tstub_DEL1\tC\t<DEL>\tDEL\tDEL\t15802\t33947957\t1\n' >> cohort_sv_pass2.tsv
    printf 'stub\t14\t99232731\t14q32.13\tstub_BND1\tG\tN[chr5:171319244[\tBND\tBND\t0\t99232731\t1\n' >> cohort_sv_pass2.tsv
    """
}
