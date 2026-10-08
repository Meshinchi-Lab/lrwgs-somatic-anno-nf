/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { PREPARE_REFERENCES } from '../subworkflows/local/prepare_references/main'
include { SAVANA_SAVANA as SAVANA } from "../modules/local/savana/savana/main"
include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main'


include { DEDUPLICATE_SAVANA } from '../subworkflows/local/deduplicate_savana'
include { MERGE_SV           } from '../subworkflows/local/merge_sv'
include { MINDA_ANNOTATIONS  } from '../subworkflows/local/minda_annotations'
include { ANNOTATIONS_SV     } from '../subworkflows/local/annotations_sv'
// Rendered TWICE -- see the two-pass block at the SV_REPORT_INDEX call site.
include { SV_REPORT_INDEX    } from '../modules/local/sv_report_index/main'
include { SV_REPORT_INDEX as SV_REPORT_INDEX_COHORT } from '../modules/local/sv_report_index/main'
include { COHORT_MERGE_SV    } from '../subworkflows/local/cohort_merge_sv/main'
include { FILTER_SV          } from '../subworkflows/local/filter_sv/main'
include { ANNOTATIONS_CNA   } from '../subworkflows/local/annotations_cna'
include { FILTER_CNA        } from '../subworkflows/local/filter_cna/main'
include { IGV_REPORTS_SV     } from '../subworkflows/local/igv_reports/main'
include { MERGE_SNV          } from '../subworkflows/local/merge_snv/main'
include { ANNOTATIONS_SNV    } from '../subworkflows/local/annotations_snv/main'

include { FILTER_SNV         } from '../subworkflows/local/filter_snv/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'


include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanopore_variant_calling_anno'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOPORE_VARIANT_CALLING_ANNO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    //
    // each sample into one item per variant type, then branch into three named channels.
    //  Each item is [ meta, [files] ] tuple convention. variant_type is embedded in meta for downstream tracing.
    //

    ch_samplesheet
        .flatMap { sample ->
            [
                [ sample.meta + [variant_type: 'sv'],  sample.svs  ],
                [ sample.meta + [variant_type: 'cnv'], sample.cnvs ],
                [ sample.meta + [variant_type: 'snv'], sample.snvs ]
            ]
        }
        .branch {
            svs:  it[0].variant_type == 'sv'
            cnvs: it[0].variant_type == 'cnv'
            snvs: it[0].variant_type == 'snv'
        }
        .set { variant_ch }

    //
    // Genome References
    //

    // fasta channel
    channel.value(file(params.fasta, checkIfExists: true))
        .map { it -> 
        def id = "${it.simpleName}"
            [  [ "id": id ] , it, [] ]
        }
        .set { fasta_ch }
    
    // boolean: whether to create a genome sizes txt file
    channel.value(params.get_sizes)
        .map { it -> it[0].toBoolean() }
        .set { sizes_out }  
    
     // boolean: whether to run the CNA annotations subworkflows 
    def run_cna_anno = params.run_cna_anno?.toString()?.toBoolean()

    // boolean: whether to run the SNV annotations subworkflows 
    def run_snv_anno = params.run_snv_anno?.toString()?.toBoolean()
   

    // index fasta file 
    SAMTOOLS_FAIDX(fasta_ch, sizes_out)

    // fasta and fai channel 
    fasta_ch
        .join(SAMTOOLS_FAIDX.out.fai , by: 0)
        .map { it -> [ it[0], it[1], it[3] ] }
        .first() 
        .set { fasta_fai } 
    
    
    //
    // SAVANA
    //

    // contig names from file
    channel.value(file(params.contigs_file, checkIfExists: true))
         .set { contigs_ch }
    
    // Branch rows on `meta.has_savana` (set upstream when a `savana_sv` column
    // is present + non-empty). Precomputed rows skip the SAVANA process and
    // feed the user-supplied VCF straight into DEDUPLICATE_SAVANA. Live rows
    // run SAVANA as before. Both paths mix back together downstream so
    // MERGE_SV/MINDA see a single unified [ meta, savana_vcf ] channel.
    variant_ch.svs
         .branch {
             precomputed: it[0].has_savana
             live:        !it[0].has_savana
         }
         .set { savana_branch }

    // savana params (live branch only — rows lacking a pre-computed VCF)
    savana_branch.live
         .map { it ->
            def int support = params.min_support
            def int a_min_reads = params.allele_min_reads
            def int cn_step_ch = params.main_cn_step_change
            def pon = params.pon_1kg

            [ it[0], it[1][0], it[1][1], pon, support, a_min_reads, cn_step_ch ]
         }
         .set { savana_ch }

    SAVANA(savana_ch, fasta_fai, contigs_ch)

    // Mix live SAVANA output with the pre-computed VCFs. svs[3] holds the
    // pre-computed savana file for `precomputed` rows (placed there by
    // PIPELINE_INITIALISATION's ch_samplesheet map).
    savana_branch.precomputed
         .map { meta, files -> [ meta, files[3] ] }
         .set { savana_precomputed_ch }

    def savana_vcf_all_ch = SAVANA.out.vcf.mix( savana_precomputed_ch )

    //
    // Dedupe SAVANA same-position INS before MINDA sees the data. Runs on
    // both live SAVANA output and pre-computed user-supplied VCFs (defensive:
    // upstream files may or may not already be deduped).
    // AnnotSV+variantconvert downstream cannot represent these (issue #20).
    //
    DEDUPLICATE_SAVANA( savana_vcf_all_ch )

    //
    // Merge the SVs
    //

    variant_ch.svs
         .map { it ->
             [ it[0], it[1][2] ]
         }
         .join( DEDUPLICATE_SAVANA.out.vcf, by: 0 )
         .map { meta, severus_vcf, savana_vcf ->
             [ meta, [ severus_vcf, savana_vcf ] ]
         }
         .set { mergesv_input_ch }
    
    channel.value(params.sv_merge_size)
         .set { sv_merge_size }

    MERGE_SV(mergesv_input_ch, sv_merge_size)

    // The canonical `candidate_genes` value channel is built after
    // ANNOTSV_SETUPUSERANNO (below) so both the file staging and its
    // downstream tuple wrapping have a single source of truth. See
    // `ch_candidate_genes_file` block after the SETUPUSERANNO invocation.

    MINDA_ANNOTATIONS(
        MERGE_SV.out.merged_sv,
        MERGE_SV.out.severus_vcf,
        MERGE_SV.out.savana_vcf
    )

    // NOTE: sv_vcf and sv_clinvar must be provided as bgzipped (.vcf.gz) + tabix-indexed (.vcf.gz.tbi) files
    def sv_custom_files = []
    if ( params.sv_vcf ) {
        sv_custom_files << file(params.sv_vcf,          checkIfExists: true)
        sv_custom_files << file("${params.sv_vcf}.tbi", checkIfExists: true)
    }
    if ( params.sv_clinvar ) {
        sv_custom_files << file(params.sv_clinvar,          checkIfExists: true)
        sv_custom_files << file("${params.sv_clinvar}.tbi", checkIfExists: true)
    }
    def sv_vcf_ch = Channel.value( sv_custom_files )

    // VEP extra plugin/annotation files (merged with any sv_vcf / sv_clinvar custom files)
    def vep_extra_files = (params.vep_extra_files ? params.vep_extra_files.collect { file(it, checkIfExists: true) } : [] ) 

    fasta_ch
        .map { meta, fa, _sz -> [ meta, fa ] }
        .set { fasta_ch2 }

    //
    // ── PREPARE_REFERENCES — stage every annotation reference once ──────
    // VEP cache, AnnotSV bundle + user BEDs (ST17 SV + ST19 CNA), the
    // candidate-gene panel, and the SNV references (ClinVar, plus CIViC
    // lifted GRCh37 -> GRCh38). See the subworkflow header for why CIViC
    // needs a liftover at all.
    //
    PREPARE_REFERENCES( fasta_ch2, params.run_snv_anno?.toString()?.toBoolean() )

    def ch_vep_cache            = PREPARE_REFERENCES.out.vep_cache
    def ch_annotsv_ann          = PREPARE_REFERENCES.out.annotsv_annotations
    def ch_candidate_genes_file = PREPARE_REFERENCES.out.candidate_genes_file
    def candidate_genes         = PREPARE_REFERENCES.out.candidate_genes

    ANNOTATIONS_SV(
        MINDA_ANNOTATIONS.out.annotated_severus_vcf,
        MINDA_ANNOTATIONS.out.savana_only_vcf,
        candidate_genes,
        params.knot_config ? file(params.knot_config, checkIfExists: true) : [],
        fasta_ch2,
        ch_vep_cache,
        sv_vcf_ch,
        vep_extra_files,
        ch_annotsv_ann
    )

    FILTER_SV(
        ANNOTATIONS_SV.out.annotsv_vcf,
        ANNOTATIONS_SV.out.vep_vcf
    )

    //
    // CNA branch: Wakhan CNV VCF to AnnotSV to KnotAnnotSV to 2-pass vembrane filter.
    // variant_ch.cnvs items are [meta, [wakhan_vcf]]; flatten to [meta, vcf].
    //
    if ( run_cna_anno ) {
        variant_ch.cnvs
            .map { meta, files -> [ meta, files[0] ] }
            .set { wakhan_vcf_ch }

        ANNOTATIONS_CNA(
            wakhan_vcf_ch,
            candidate_genes,
            params.knot_config ? file(params.knot_config, checkIfExists: true) : [],
            ch_annotsv_ann
        )

        FILTER_CNA( ANNOTATIONS_CNA.out.annotsv_vcf, fasta_ch2 )
    }


    //
    // ── SNV intersection (ClairS-TO ∩ DeepSomatic) ─────────────────────
    // variant_ch.snvs channel are [ meta(+variant_type:'snv'), [clairsto_vcf, deepsomatic_vcf] ] and converted to [ meta_without_variant_type, clairsto_vcf, deepsomatic_vcf ] 
    //
    if ( run_snv_anno ) {
        snv_input_ch = variant_ch.snvs
            .map { meta, files ->
                def m = meta.findAll { it.key != 'variant_type' }
                [ m, files[0], files[1] ]
            }
        MERGE_SNV( snv_input_ch, fasta_fai, PREPARE_REFERENCES.out.rescue_bed )


        //
        // ── ANNOTATIONS_SNV — 
        //
        ANNOTATIONS_SNV(
            MERGE_SNV.out.vcf,
            MERGE_SNV.out.tbi,
            PREPARE_REFERENCES.out.clinvar_vcf,
            PREPARE_REFERENCES.out.clinvar_tbi,
            PREPARE_REFERENCES.out.civic_vcf,
            PREPARE_REFERENCES.out.civic_tbi,
            ch_vep_cache,
            fasta_ch2
        )

        //
        // ── FILTER_SNV —
        //
        FILTER_SNV(
            ANNOTATIONS_SNV.out.vcf,
            ANNOTATIONS_SNV.out.tbi,
            // `params.candidate_genes` (undeclared) was used here until 2026-09-17. It
            // resolved to null, so FILTER_SNV's clause D dropped its gene predicate
            // entirely and admitted every gene instead of the candidate panel. The
            // canonical source is PREPARE_REFERENCES.out.candidate_genes_file, defined
            // above as `ch_candidate_genes_file` (raw path, or [] when unset).
            ch_candidate_genes_file
        )
    }

    // Per-sample IGV interactive HTML scoped to the pass2 prioritized SVs.
    // variant_ch.svs items are [ meta, [bam, bai, severus_vcf, ...] ]; take
    // the first two elements as the BAM+BAI pair
    variant_ch.svs
        .map { meta, files -> [ meta, files[0], files[1] ] }
        .set { sample_bam_ch }

    // Guard channels: substitute Channel.empty() when the CNA or SNV branch
    // is disabled. CNA channels are no longer consumed by IGV_REPORTS_SV 
    def cna_tsv_pass2_ch = run_cna_anno ? FILTER_CNA.out.tsv_pass2 : Channel.empty()
    def snv_tsv_ch       = run_snv_anno ? FILTER_SNV.out.tsv       : Channel.empty()

    // Data dictionaries emitted by CREATE_DATA_DICT_{SV,CNA,SNV} inside each
    // filter subworkflow. Guarded the same way as the TSV channels so a
    // disabled branch collapses to an empty channel (SV_REPORT_INDEX's stage
    // spec is optional-list-friendly via ifEmpty([]) below).
    def sv_dict_ch  = FILTER_SV.out.dict
    def cna_dict_ch = run_cna_anno ? FILTER_CNA.out.dict : Channel.empty()
    def snv_dict_ch = run_snv_anno ? FILTER_SNV.out.dict : Channel.empty()

    // IGV reports consume the SV pass2 VCF, and (when enabled) the FILTER_SNV VCF for merging into the SMALL bucket.
    def snv_vcf_igv_ch = run_snv_anno ? FILTER_SNV.out.vcf : Channel.empty()
    def snv_tbi_igv_ch = run_snv_anno ? FILTER_SNV.out.tbi : Channel.empty()

    IGV_REPORTS_SV(
        sample_bam_ch,
        FILTER_SV.out.vcf_pass2,
        FILTER_SV.out.tbi_pass2,
        snv_vcf_igv_ch,
        snv_tbi_igv_ch,
        fasta_fai,
        ch_candidate_genes_file
    )

    //
    // ── SV_REPORT_INDEX — TWO-PASS RENDER ───────────────────────────────
    // driver_tier exists only in the Quarto layer, so COHORT_MERGE_SV consumes
    // this process's output while its own results belong in the same report — a
    // genuine dependency cycle for a single render. Rather than lift ~300 lines
    // of interdependent tiering logic out of the qmd into a module, the SAME qmd
    // is rendered twice:
    //
    //   pass 1  SV_REPORT_INDEX        → cohort_sv_pass2.tsv + "pending" cohort card
    //   then    COHORT_MERGE_SV        → cohort_sv_recurrent.tsv + merged VCF
    //   pass 2  SV_REPORT_INDEX_COHORT → same report, cohort section populated
    //
    // The tiering code never moves between passes, so the tiers pass 2 renders
    // ARE the tiers the merge consumed — identical by construction, not by
    // assumption. Both passes publish to results/sv_report/; pass 2 necessarily
    // finishes later (it depends on pass 1 through COHORT_MERGE_SV), so its
    // index.html wins. With run_cohort_merge_sv=false there is no pass 2 and
    // pass 1's report — cohort card marked unavailable — is the published one.
    //
    // The inputs are bound to locals so both passes are fed the SAME channels;
    // DSL2 forks a channel across multiple consumers (ch_candidate_genes_file
    // already feeds three).
    //
    def rpt_annotsv_tsv = ANNOTATIONS_SV.out.tsv.map { _meta, tsv -> tsv }.collect()
    def rpt_knot_html   = ANNOTATIONS_SV.out.html.map { _meta, html -> html }.collect()
    def rpt_knot_xl     = ANNOTATIONS_SV.out.xl.map { _meta, xlsm -> xlsm }.collect().ifEmpty( [] )
    def rpt_vep_html    = ANNOTATIONS_SV.out.vep_report.map { _meta, html -> html }.collect().ifEmpty( [] )
    def rpt_igv_html    = IGV_REPORTS_SV.out.report.map { _meta, html -> html }.collect().ifEmpty( [] )
    def rpt_sv_pass2    = FILTER_SV.out.tsv_pass2.map { _meta, tsv -> tsv }.collect().ifEmpty( [] )
    def rpt_cna_pass2   = cna_tsv_pass2_ch.map { _meta, tsv -> tsv }.collect().ifEmpty( [] )
    def rpt_snv_tsv     = snv_tsv_ch.map { _meta, tsv -> tsv }.collect().ifEmpty( [] )
    def rpt_metadata    = params.sample_metadata ? file(params.sample_metadata, checkIfExists: true) : []
    def rpt_ensdb       = params.ensembl_db      ? file(params.ensembl_db, checkIfExists: true)      : []
    // AnnotSV annotations dir — symlinked into the SV_REPORT_INDEX work dir so
    // the qmd can glob `Annotations_Human/Users/**/*.header.tsv` to derive the ST
    // regex used by parse_annotsv(). Extract just the dir path (drop the meta)
    // from the ch_annotsv_ann value-channel.
    def rpt_annotsv_ann = ch_annotsv_ann.map { _meta, dir -> dir }
    // Per-sample VCF header dictionaries → data_dict_{sv,cna,snv}/ subdirs.
    // Union'd inside the qmd across samples (same field_name → same description,
    // so dedup keeps a single row per definition).
    def rpt_dict_sv     = sv_dict_ch .map { _meta, tsv -> tsv }.collect().ifEmpty( [] )
    def rpt_dict_cna    = cna_dict_ch.map { _meta, tsv -> tsv }.collect().ifEmpty( [] )
    def rpt_dict_snv    = snv_dict_ch.map { _meta, tsv -> tsv }.collect().ifEmpty( [] )

    SV_REPORT_INDEX(
        file("${projectDir}/bin/sv_report_index.qmd"),
        rpt_annotsv_tsv,
        rpt_knot_html,
        rpt_knot_xl,
        rpt_vep_html,
        rpt_igv_html,
        rpt_sv_pass2,
        rpt_cna_pass2,
        rpt_snv_tsv,
        rpt_metadata,
        // Canonical candidate-gene panel — sourced from SETUPUSERANNO.out.genes
        // (or the workflow-level fallback when SETUPUSERANNO is skipped).
        ch_candidate_genes_file,
        rpt_ensdb,
        rpt_annotsv_ann,
        rpt_dict_sv,
        rpt_dict_cna,
        rpt_dict_snv,
        [],                     // cohort_recurrent_tsv — pass 1: merge has not run
        [],                     // cohort_merged_vcf    — pass 1: merge has not run
        // Pass 1's report is overwritten by pass 2 whenever the cohort merge
        // runs, so skip the per-sample DISCO rendering — the most expensive part
        // of the render — and let pass 2 produce it. With the merge disabled
        // there is no pass 2, so pass 1 must render DISCO normally.
        params.run_cohort_merge_sv?.toString()?.toBoolean() ? false : params.disco_enable
    )

    //
    // ── COHORT_MERGE_SV — recurrent SVs across the cohort ───────────────
    // Runs AFTER SV_REPORT_INDEX because driver_tier is computed in the Quarto layer
    // and exists in no VCF. The tiered TSV supplies the REGIONS; the records searched
    // are the PRE-FILTER VEP-annotated VCFs, so an SV called confidently in one sample
    // and filtered out in another is rescued rather than lost.
    //
    // NOTE: SV_REPORT_INDEX.out.tsv_sv is a bare `path` (and `optional: true`), so it is
    // wrapped with a cohort meta here to match COLLECT_COORDINATES' tuple input.
    // ANNOTATIONS_SV.out.vep_vcf is deliberately used instead of FILTER_SV's pass2 VCFs:
    // restricting pass2 to pass2 would make rescue impossible.
    //
    if ( params.run_cohort_merge_sv?.toString()?.toBoolean() ) {
        COHORT_MERGE_SV(
            SV_REPORT_INDEX.out.tsv_sv.map { tsv -> [ [ id: 'cohort' ], tsv ] },
            ANNOTATIONS_SV.out.vep_vcf.join( ANNOTATIONS_SV.out.vep_tbi, by: 0 ),
            fasta_ch2,
            fasta_fai.map { meta, _fa, fai -> [ meta, fai ] }
        )

        //
        // Pass 2 — the SAME qmd, re-rendered with the cohort results staged.
        // COHORT_MERGE_SV skips Jasmine entirely below N=2 (its N=1 guard), so
        // both emits can legitimately be empty. .ifEmpty( [] ) keeps pass 2
        // scheduled in that case and the qmd falls back to the "unavailable"
        // card, instead of the whole render silently dropping out of the DAG.
        //
        SV_REPORT_INDEX_COHORT(
            file("${projectDir}/bin/sv_report_index.qmd"),
            rpt_annotsv_tsv,
            rpt_knot_html,
            rpt_knot_xl,
            rpt_vep_html,
            rpt_igv_html,
            rpt_sv_pass2,
            rpt_cna_pass2,
            rpt_snv_tsv,
            rpt_metadata,
            ch_candidate_genes_file,
            rpt_ensdb,
            rpt_annotsv_ann,
            rpt_dict_sv,
            rpt_dict_cna,
            rpt_dict_snv,
            COHORT_MERGE_SV.out.recurrent .map { _meta, tsv -> tsv }.ifEmpty( [] ),
            COHORT_MERGE_SV.out.merged_vcf.map { _meta, vcf -> vcf }.ifEmpty( [] ),
            params.disco_enable     // pass 2 is the published report: render DISCO
        )
    }

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    /*
    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'lrwgs-somatic-anno-nf_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'lrwgs-somatic-anno-nf'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    */
    
    emit: 
    // multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
