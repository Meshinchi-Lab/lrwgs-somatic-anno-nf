//
// Annotate merged SVs with AnnotSV and produce a KnotAnnotSV interactive report
//

include { BCFTOOLS_CONCAT                          } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_MERGED  } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_SORT  as BCFTOOLS_SORT_MERGED   } from '../../../modules/nf-core/bcftools/sort/main'

include { ANNOTSV_ANNOTSV                          } from '../../../modules/local/annotsv/annotsv/main'
include { BCFTOOLS_REHEADER_SAVANA_FMT             } from '../../../modules/local/bcftools/reheader_savana_fmt/main'
include { BCFTOOLS_CLEAN_ANNOTSV as BCFTOOLS_CLEAN_ANNOTSV_SV } from '../../../modules/local/bcftools/clean_annotsv/main'
include { ENSEMBLVEP_VEP                           } from '../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_FILTERVEP                     } from '../../../modules/nf-core/ensemblvep/filtervep/main'
include { SPLIT_VEP_TO_INFO                        } from '../../../modules/local/split_vep_to_info/main'
// include { SVDB_QUERY                            } from '../../../modules/nf-core/svdb/query/main'
include { VEMBRANE_TABLE as VEMBRANE_TABLE_CALLER  } from '../../../modules/nf-core/vembrane/table/main'
include { ADD_CALLER                               } from '../../../modules/local/add_caller/main'
include { KNOTANNOTSV as KNOTANNOTSV_HTML          } from '../../../modules/nf-core/knotannotsv/main'
include { KNOTANNOTSV as KNOTANNOTSV_XL            } from '../../../modules/nf-core/knotannotsv/main'

workflow ANNOTATIONS_SV {

    take:
    annotated_severus_ch  // channel: [ val(meta), path(vcf.gz) ] — Severus VCF with CALLER/SAVANA INFO
    savana_only_vcf_ch    // channel: [ val(meta), path(vcf.gz) ] — SAVANA-only calls, sorted
    candidate_genes       // val: path
    knot_config           // val: path or [] — optional custom config_AnnotSV.yaml
    ch_fasta              // channel: [ val(meta), path(fasta) ] — reference FASTA for VEP
    ch_vep_cache          // channel: [ val(meta2), path(cache) ] — VEP cache dir ([] when not used)
    sv_vcf_ch             // Channel.value( list<path> ) — [sv_vcf, sv_vcf_tbi, clinvar, clinvar_tbi] or [] when unused
    vep_extra_files       // val: list<path> — extra files staged for VEP plugins ([] if none)
    ch_annotsv_ann        // channel: [ val(meta), path(annotations_dir) ] — staged by main wf (shared with ANNOTATIONS_CNA)

    main:

    // Step 2: concatenate annotated Severus + SAVANA-only VCFs
    // tbi=[] lets the module run tabix internally for any unindexed input
    // Step 2: optionally concatenate annotated Severus + SAVANA-only VCFs.
    // When params.add_savana_sv is true, BCFTOOLS_CONCAT joins both callers'
    // calls before AnnotSV. When false (default), only the annotated Severus
    // VCF is carried forward — this keeps the downstream KnotAnnotSV HTML
    // small enough to open in a browser by dropping the SAVANA-only records.
    def ch_pre_sort
    if ( params.add_savana_sv ) {
        annotated_severus_ch
            .join( savana_only_vcf_ch, by: 0 )
            .map { meta, sev_ann, sav_only -> [ meta, [ sev_ann, sav_only ], [] ] }
            .set { ch_concat_input }

        BCFTOOLS_CONCAT( ch_concat_input )
        ch_pre_sort = BCFTOOLS_CONCAT.out.vcf
    } else {
        ch_pre_sort = annotated_severus_ch
    }

    // Include another optional subworkflow, so that the full cohort of samples are merged into a single VCF
    // then finish the same process of annotations and filtering?

    // Step 3: sort and index the chosen VCF (Severus-only or Severus+SAVANA concat)
    BCFTOOLS_SORT_MERGED(  ch_pre_sort )
    BCFTOOLS_INDEX_MERGED( BCFTOOLS_SORT_MERGED.out.vcf )

    // Step 4: AnnotSV annotation
    // ANNOTSV_ANNOTSV input[0]: tuple val(meta), path(sv_vcf), path(sv_vcf_index), path(candidate_small_variants)
    BCFTOOLS_SORT_MERGED.out.vcf
        .join( BCFTOOLS_INDEX_MERGED.out.tbi, by: 0 )
        .map { meta, vcf, tbi -> [ meta, vcf, tbi, [] ] }
        .set { ch_annotsv_input }

    ANNOTSV_ANNOTSV(
        ch_annotsv_input,
        ch_annotsv_ann,
        candidate_genes,
        Channel.value( [ [:], [] ] ),   // no false_positive_snv
        Channel.value( [ [:], [] ] )    // no gene_transcripts
    )

    // Step 5: SVDB population frequency query
    // SVDB_QUERY(
    //     ANNOTSV_ANNOTSV.out.vcf,
    //     [],           // in_occs
    //     [],           // in_frqs
    //     [],           // out_occs
    //     [],           // out_frqs
    //     svdb_vcf_dbs,
    //     []            // no bedpe_dbs
    // )

    // Step 5b: inject missing SAVANA FORMAT header definitions (VAF, hVAF, DR, DV).
    // AnnotSV does not forward FORMAT ##FORMAT lines from its input VCF, so these
    // fields are present in records but undefined in the header, causing pysam /
    // vembrane to fail.  bcftools annotate --header-lines is additive (appends only).
    BCFTOOLS_REHEADER_SAVANA_FMT( ANNOTSV_ANNOTSV.out.vcf )

    // Step 5c: normalise ACMG_class + AnnotSV_ranking_score to clean scalars
    // (Number=1,Type=Integer / Type=Float). AnnotSV emits ACMG_class as a
    // compound `<split>,full=<full>` string with Number=. header, forcing
    // downstream vembrane / KnotAnnotSV / the R report to string-parse it.
    // Cleaning here (before FILTER_SV) means every consumer sees typed scalars.
    BCFTOOLS_CLEAN_ANNOTSV_SV( BCFTOOLS_REHEADER_SAVANA_FMT.out.vcf )

    // Step 6a: extract CALLER (and any other small INFO subset) from the merged VCF
    // via vembrane/table. Runs in parallel with ANNOTSV_ANNOTSV — same input VCF,
    // no AnnotSV dependency. Produces a 2-column TSV (ID + CALLER) used downstream
    // to replace AnnotSV's bloated full-INFO column in the KnotAnnotSV HTML.
    // VAF, hVAF, DR, DV come from the per-sample FORMAT column from Severus 
    VEMBRANE_TABLE_CALLER(
        BCFTOOLS_SORT_MERGED.out.vcf,
        'ID, INFO.get("CALLER", ""), FORMAT["VAF"][SAMPLE] if "VAF" in FORMAT else "", FORMAT["hVAF"][SAMPLE] if "hVAF" in FORMAT else "", FORMAT["DR"][SAMPLE] if "DR" in FORMAT else "", FORMAT["DV"][SAMPLE] if "DV" in FORMAT else ""'
    )

    // Step 6b: merge AnnotSV TSV with the caller lookup on ID; drop the full INFO
    // column so KnotAnnotSV emits a slim HTML/Excel report. ADD_CALLER also
    // computes a `Candidate_genes` column by intersecting each row's
    // `Gene_name` with the user's candidate-gene list (the same file passed
    // to AnnotSV via `-candidateGenesFile`, which only affects ranking and
    // does not emit a per-row hit list of its own).
    ANNOTSV_ANNOTSV.out.tsv
        .join( VEMBRANE_TABLE_CALLER.out.table, by: 0 )
        .set { ch_add_caller_input }

    // Extract the raw path from the `candidate_genes` take input (tuple
    // form `[ meta, path_or_empty ]`) so ADD_CALLER can stage it and pass
    // --candidate-genes to add_caller.py. Sourcing from the take input
    // (which flows from ANNOTSV_SETUPUSERANNO.out.genes at the workflow
    // level) instead of re-reading params.candidateGenesFile here keeps
    // the panel file's provenance a single-DAG-source.
    def candidate_genes_file = candidate_genes.map { _meta, f -> f ?: [] }

    ADD_CALLER( ch_add_caller_input, candidate_genes_file )

    // Step 6c: KnotAnnotSV — emit BOTH the interactive HTML and the Excel
    // (.xlsm) report from the same slim TSV. The upstream nf-core module
    // produces one or the other per invocation (the boolean knot_out_xl
    // toggles between knotAnnotSV.pl and knotAnnotSV2XL.pl)
    ADD_CALLER.out.tsv_html
        .map { meta, tsv -> [ meta, tsv, false, knot_config ] }
        .set { ch_knot_input_html }

    ADD_CALLER.out.tsv
        .map { meta, tsv -> [ meta, tsv, true,  knot_config ] }
        .set { ch_knot_input_xl }

    KNOTANNOTSV_HTML( ch_knot_input_html )
    KNOTANNOTSV_XL(   ch_knot_input_xl   )

    // Step 7: VEP annotation of the AnnotSV VCF
    // ENSEMBLVEP_VEP input[0]: [ meta, vcf, custom_extra_files ]
    // sv_vcf_ch is Channel.value(list) — always emits one item (possibly []) and .combine() always produces output and the list drops cleanly into
    // VEP's path(custom_extra_files), stages all files into the work dir.
    def ch_vep_vcf    = channel.empty()
    def ch_vep_tbi    = channel.empty()
    def ch_vep_report = channel.empty()
    def run_vep = params.sv_anno_vep ? params.sv_anno_vep.toString().toBoolean() : false

    // sv_anno_vep accepts either boolean or string ('true'/'false') in the params block.
    if ( run_vep ) {
        // .combine() spreads Channel.value(list) into individual tuple elements,
        // so .drop(2) re-collects everything after [meta, vcf] back into a list.
        // VEP consumes the ACMG-cleaned VCF so its output preserves the typed
        // Number=1 ACMG_class + AnnotSV_ranking_score alongside CSQ additions.
        BCFTOOLS_CLEAN_ANNOTSV_SV.out.vcf
            .combine( sv_vcf_ch )
            .map { it -> [ it[0], it[1], it.drop(2) ] }
            .set { ch_vep_input }
        ENSEMBLVEP_VEP(
            ch_vep_input,
            params.genome,
            params.species,
            params.cache_version,
            ch_vep_cache,
            ch_fasta,
            vep_extra_files
        )
        // Step 7b: extract named CSQ subfields (gnomAD_SV, gnomAD_SV_AF by default)into INFO feilds
        SPLIT_VEP_TO_INFO( ENSEMBLVEP_VEP.out.vcf )
        ch_vep_vcf    = SPLIT_VEP_TO_INFO.out.vcf
        ch_vep_tbi    = SPLIT_VEP_TO_INFO.out.tbi
        // ENSEMBLVEP_VEP.out.report shape: [meta, task.process, 'ensemblvep', path(*.html)]
        ch_vep_report = ENSEMBLVEP_VEP.out.report
                            .map { meta, _proc, _tool, html -> [ meta, html ] }

        // // Step 7b: convert VEP VCF to per-SV tab summary (--pick to one row per variant)
        // ENSEMBLVEP_FILTERVEP(
        //     ch_vep_vcf,
        //     [],     // no feature_file
        //     'tab'   // output extension
        // )
        // ch_vep_tab = ENSEMBLVEP_FILTERVEP.out.output
    }

    emit:
    // vcf          = SVDB_QUERY.out.vcf        // final annotated VCF per sample
    tsv          = ANNOTSV_ANNOTSV.out.tsv            // AnnotSV full TSV
    annotsv_vcf  = BCFTOOLS_CLEAN_ANNOTSV_SV.out.vcf  // [ meta, *.acmg_clean.vcf.gz ] — SAVANA FORMAT headers restored + ACMG_class/AnnotSV_ranking_score normalised to typed scalars
    html         = KNOTANNOTSV_HTML.out.html          // [ meta, *.html  ] — KnotAnnotSV interactive HTML report
    xl           = KNOTANNOTSV_XL.out.xl              // [ meta, *.xlsm  ] — KnotAnnotSV Excel report
    vep_vcf      = ch_vep_vcf                         // [ meta, vcf.gz ] — VEP annotated VCF (empty when VEP skipped)
    vep_tbi      = ch_vep_tbi                         // [ meta, vcf.gz.tbi ] — index for the above (empty when VEP skipped)
    vep_report   = ch_vep_report                      // [ meta, path(*.html) ] — VEP summary HTML (empty when VEP skipped)
    versions     = channel.empty()
}
