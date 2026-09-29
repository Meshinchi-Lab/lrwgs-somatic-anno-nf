// subworkflows/local/annotations_cna/main.nf
//
// Annotate Wakhan CNV VCF with AnnotSV and produce KnotAnnotSV HTML+XL reports.
// Mirrors ANNOTATIONS_SV. Adds BCFTOOLS_ADD_CALLER_TAG upstream of sort to include
// a constant `CALLER=WAKHAN` INFO tag on every record for consistency and if a 2nd caller is included later. 
//

include { BCFTOOLS_REHEADER                        } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_RENAME_IDS                      } from '../../../modules/local/bcftools/rename_ids/main'
include { BCFTOOLS_ADD_CALLER_TAG                  } from '../../../modules/local/bcftools/add_caller_tag/main'
include { BCFTOOLS_SORT  as BCFTOOLS_SORT_CNA      } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CNA     } from '../../../modules/nf-core/bcftools/index/main'

include { ANNOTSV_ANNOTSV                          } from '../../../modules/local/annotsv/annotsv/main'
include { BCFTOOLS_REHEADER_WAKHAN_FMT             } from '../../../modules/local/bcftools/reheader_wakhan_fmt/main'
include { BCFTOOLS_CLEAN_ANNOTSV as BCFTOOLS_CLEAN_ANNOTSV_CNA } from '../../../modules/local/bcftools/clean_annotsv/main'

include { VEMBRANE_TABLE as VEMBRANE_TABLE_CALLER  } from '../../../modules/nf-core/vembrane/table/main'
include { ADD_CALLER                               } from '../../../modules/local/add_caller/main'
include { KNOTANNOTSV as KNOTANNOTSV_HTML          } from '../../../modules/nf-core/knotannotsv/main'
include { KNOTANNOTSV as KNOTANNOTSV_XL            } from '../../../modules/nf-core/knotannotsv/main'

workflow ANNOTATIONS_CNA {

    take:
    wakhan_vcf_ch    // channel: [ val(meta), path(vcf.gz) ] — Wakhan CN-integers VCF
    candidate_genes  // val: path or [] — consumed by ADD_CALLER --candidate-genes
    knot_config      // val: path or [] — optional custom config_AnnotSV.yaml
    ch_annotsv_ann   // channel.value([ [:], path(annotations_dir) ]) — staged in main wf

    main:

    // Step 0: rewrite VCF IDs to coordinate-free sequential strings (wakhan_1,wakhan_2, ...) and change the default column name "Sample" to the sample_id from the sample_sheet input
    // AnnotSV's variantconvert shifts POS by +1 and rewrites coordinates embedded in IDs, which causes down-stream merge error;
    //  sequential IDs survive unchanged so them downstream ADD_CALLER pandas merge on ID.
            // Rename the SEVERUS sample column to meta.id.
    // bcftools reheader --samples takes a plain-text file with the new sample name(s).
    wakhan_vcf_ch
        .map { meta, vcf ->
            def sf = File.createTempFile("${meta.id}.${meta.variant_type}.samples", ".txt", new File("${workflow.workDir}"))
            sf.deleteOnExit()
            sf.text = "${meta.id}\n"
            [ meta, vcf, [], file(sf) ]
        }
        .set { ch_reheader_input }

    BCFTOOLS_REHEADER( ch_reheader_input, Channel.value( [ [:], [] ] ) )

    
    def id_prefix = 'wakhan'
    BCFTOOLS_RENAME_IDS( BCFTOOLS_REHEADER.out.vcf, id_prefix )



    // Step 1: tag every record with CALLER=WAKHAN
    add_caller_input_ch = BCFTOOLS_RENAME_IDS.out.vcf
        .join( BCFTOOLS_RENAME_IDS.out.tbi, by: 0 )
    BCFTOOLS_ADD_CALLER_TAG( add_caller_input_ch, 'WAKHAN', id_prefix)

    // Step 2: sort + index the tagged VCF
    BCFTOOLS_SORT_CNA(  BCFTOOLS_ADD_CALLER_TAG.out.vcf.map { meta, vcf, tbi -> [ meta, vcf ] } )
    BCFTOOLS_INDEX_CNA( BCFTOOLS_SORT_CNA.out.vcf )

    // Step 3: AnnotSV annotation
    // ANNOTSV_ANNOTSV input[0]: tuple val(meta), path(sv_vcf), path(sv_vcf_index), path(candidate_small_variants)
    BCFTOOLS_SORT_CNA.out.vcf
        .join( BCFTOOLS_INDEX_CNA.out.tbi, by: 0 )
        .map { meta, vcf, tbi -> [ meta, vcf, tbi, [] ] }
        .set { ch_annotsv_input }

    ANNOTSV_ANNOTSV(
        ch_annotsv_input,
        ch_annotsv_ann,
        candidate_genes,
        Channel.value( [ [:], [] ] ),   // no false_positive_snv
        Channel.value( [ [:], [] ] )    // no gene_transcripts
    )

    // Step 4: restore Wakhan FORMAT headers (AnnotSV/variantconvert strips them).
    // The original Wakhan VCF is the source for FORMAT header content.
    ANNOTSV_ANNOTSV.out.vcf
        .join( wakhan_vcf_ch, by: 0 )
        .set { ch_reheader_input }

    // id_prefix is the shared Wakhan record-ID prefix set upstream by
    // BCFTOOLS_RENAME_IDS. REHEADER_WAKHAN_FMT needs it to scope the
    // 0↔1 POS-shift sed (workaround for bcftools annotate refusing POS=0
    // — see module header comment) to Wakhan records only.
    BCFTOOLS_REHEADER_WAKHAN_FMT( ch_reheader_input, id_prefix )

    // Normalise ACMG_class + AnnotSV_ranking_score to typed scalars before
    // FILTER_CNA so vembrane, KnotAnnotSV, and the R report all see clean
    // integers/floats instead of AnnotSV's `<split>,full=<full>` compound
    // encoding. Aliased alongside the SV branch (BCFTOOLS_CLEAN_ANNOTSV_SV)
    // so per-subworkflow publishDir catch-alls place the outputs correctly.
    BCFTOOLS_CLEAN_ANNOTSV_CNA( BCFTOOLS_REHEADER_WAKHAN_FMT.out.vcf )

    // Step 5: caller-lookup TSV (CALLER + Wakhan FORMAT fields)
    // Reads the sorted+tagged VCF (before AnnotSV), where the FORMAT headers are
    // still defined. AnnotSV have strips them.
    VEMBRANE_TABLE_CALLER(
        BCFTOOLS_SORT_CNA.out.vcf,
        'ID, INFO.get("CALLER", ""), FORMAT["TCN"][SAMPLE] if "TCN" in FORMAT else "", FORMAT["CN1"][SAMPLE] if "CN1" in FORMAT else "", FORMAT["CN2"][SAMPLE] if "CN2" in FORMAT else "", FORMAT["CNQ1"][SAMPLE] if "CNQ1" in FORMAT else "", FORMAT["CNQ2"][SAMPLE] if "CNQ2" in FORMAT else "", FORMAT["COV1"][SAMPLE] if "COV1" in FORMAT else "", FORMAT["COV2"][SAMPLE] if "COV2" in FORMAT else ""'
    )

    // Step 6: join AnnotSV TSV with caller-lookup; drop the bloated full-INFO
    // column; compute Candidate_genes by intersecting Gene_name with the candidate-genes list
    ANNOTSV_ANNOTSV.out.tsv
        .join( VEMBRANE_TABLE_CALLER.out.table, by: 0 )
        .set { ch_add_caller_input }

    // Extract the raw path from the `candidate_genes` take input (tuple
    // form `[ meta, path_or_empty ]`) — canonical source is
    // ANNOTSV_SETUPUSERANNO.out.genes, forwarded through the workflow. This
    // replaces the previous direct params.candidateGenesFile read so the
    // panel file's provenance stays a single-DAG-source.
    def candidate_genes_file = candidate_genes.map { _meta, f -> f ?: [] }

    ADD_CALLER( ch_add_caller_input, candidate_genes_file )

    // Step 7: KnotAnnotSV — emit BOTH the interactive HTML and the Excel (.xlsm) report from the same  TSV. 
    ADD_CALLER.out.tsv_html
        .map { meta, tsv -> [ meta, tsv, false, knot_config ] }
        .set { ch_knot_input_html }

    ADD_CALLER.out.tsv
        .map { meta, tsv -> [ meta, tsv, true, knot_config ] }
        .set { ch_knot_input_xl }

    KNOTANNOTSV_HTML( ch_knot_input_html )
    KNOTANNOTSV_XL(   ch_knot_input_xl )

    emit:
    tsv         = ANNOTSV_ANNOTSV.out.tsv                  // [ meta, *.tsv ]
    annotsv_vcf = BCFTOOLS_CLEAN_ANNOTSV_CNA.out.vcf     // [ meta, *.acmg_clean.vcf.gz ] — Wakhan FMT headers restored + ACMG_class/AnnotSV_ranking_score normalised to typed scalars; input to FILTER_CNA
    html        = KNOTANNOTSV_HTML.out.html                // [ meta, *.html ]
    xl          = KNOTANNOTSV_XL.out.xl                    // [ meta, *.xlsm ]
    versions    = channel.empty()
}
