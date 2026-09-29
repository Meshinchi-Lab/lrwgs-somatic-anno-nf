// subworkflows/local/merge_snv/main.nf
//
// Intersect ClairS-TO + DeepSomatic SNV calls (2-of-2 strict agreement using
// `bcftools isec -n+2 -c none`), then enrich the consensus VCF with per-caller
// VAF and DP as INFO tags (CLAIRSTO_VAF, CLAIRSTO_DP, DEEPSOMATIC_VAF,
// DEEPSOMATIC_DP) plus an INFO/CALLER provenance tag. Warning-only QC step
// compares the final consensus VCF to the isec ground truth.

include { BCFTOOLS_NORM         as BCFTOOLS_NORM_CLAIRSTO     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM         as BCFTOOLS_NORM_DEEPSOMATIC  } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW         as BCFTOOLS_VIEW_PASS_CLAIRSTO    } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW         as BCFTOOLS_VIEW_PASS_DEEPSOMATIC } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ISEC                                  } from '../../../modules/local/bcftools/isec/main'

include { BCFTOOLS_RENAME_CALLER_FIELDS as RENAME_CLAIRSTO    } from '../../../modules/local/bcftools/rename_caller_fields/main'
include { BCFTOOLS_RENAME_CALLER_FIELDS as RENAME_DEEPSOMATIC } from '../../../modules/local/bcftools/rename_caller_fields/main'
include { BCFTOOLS_VIEW_TO_INFO         as VIEW_TO_INFO_CLAIRSTO    } from '../../../modules/local/bcftools/view_to_info/main'
include { BCFTOOLS_VIEW_TO_INFO         as VIEW_TO_INFO_DEEPSOMATIC } from '../../../modules/local/bcftools/view_to_info/main'

include { BCFTOOLS_ANNOTATE    as BCFTOOLS_ANNOTATE_CROSS_CALLER } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ADD_CALLER_TAG                              } from '../../../modules/local/bcftools/add_caller_tag/main'
include { BCFTOOLS_SORT        as BCFTOOLS_SORT_MERGE_SNV       } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_REHEADER                                    } from '../../../modules/nf-core/bcftools/reheader/main'
include { QC_MERGE_SNV                                          } from '../../../modules/local/qc/merge_snv/main'

include { BCFTOOLS_VIEW  as RESCUE_VIEW_CLAIRSTO    } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW  as RESCUE_VIEW_DEEPSOMATIC } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_RENAME_CALLER_FIELDS as RESCUE_RENAME_CLAIRSTO    } from '../../../modules/local/bcftools/rename_caller_fields/main'
include { BCFTOOLS_RENAME_CALLER_FIELDS as RESCUE_RENAME_DEEPSOMATIC } from '../../../modules/local/bcftools/rename_caller_fields/main'
include { BCFTOOLS_VIEW_TO_INFO as RESCUE_INFO_CLAIRSTO    } from '../../../modules/local/bcftools/view_to_info/main'
include { BCFTOOLS_VIEW_TO_INFO as RESCUE_INFO_DEEPSOMATIC } from '../../../modules/local/bcftools/view_to_info/main'
include { BCFTOOLS_ADD_CALLER_TAG as RESCUE_TAG_CLAIRSTO    } from '../../../modules/local/bcftools/add_caller_tag/main'
include { BCFTOOLS_ADD_CALLER_TAG as RESCUE_TAG_DEEPSOMATIC } from '../../../modules/local/bcftools/add_caller_tag/main'
include { BCFTOOLS_ANNOTATE as RESCUE_FLAG_CLAIRSTO    } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as RESCUE_FLAG_DEEPSOMATIC } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_SUBTRACT_VCF as RESCUE_SUBTRACT_CLAIRSTO    } from '../../../modules/local/bcftools/subtract_vcf/main'
include { BCFTOOLS_SUBTRACT_VCF as RESCUE_SUBTRACT_DEEPSOMATIC } from '../../../modules/local/bcftools/subtract_vcf/main'
include { BCFTOOLS_REHEADER as RESCUE_REHEADER_CONSENSUS   } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as RESCUE_REHEADER_CLAIRSTO    } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as RESCUE_REHEADER_DEEPSOMATIC } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_CONCAT as MERGE_SNV_CONCAT_RESCUE } from '../../../modules/nf-core/bcftools/concat/main'

workflow MERGE_SNV {
    take:
    snv_ch       // channel: [ val(meta), path(clairsto_vcf), path(deepsomatic_vcf) ]
    fasta_ch     // value:   [ [id], path(fasta), path(fai) ]   — same shape as fasta_fai
    rescue_bed     // value: [ meta, bed ] or [ [:], [] ] from PREPARE_REFERENCES

    main:

    // Fan into per-caller inputs (input vcf paired with empty tbi placeholder;
    // bcftools/norm tolerates missing tbi and will recompute).
    clairsto_in_ch    = snv_ch.map { meta, c, d -> [ meta, c, [] ] }
    deepsomatic_in_ch = snv_ch.map { meta, c, d -> [ meta, d, [] ] }

    // ── Per-caller normalization (split multi-allelics, left-align, trim) ──
    // bcftools/norm second input is `tuple val(meta2), path(fasta)`; we use the
    // existing fasta_ch which already has [meta, fasta, fai] shape — drop fai.
    fasta_only_ch = fasta_ch.map { m, fa, fi -> [ m, fa ] }
    BCFTOOLS_NORM_CLAIRSTO(    clairsto_in_ch,    fasta_only_ch )
    BCFTOOLS_NORM_DEEPSOMATIC( deepsomatic_in_ch, fasta_only_ch )

    // ── PASS-filter each ────────────────────────────────────────────────
    clairsto_view_ch    = BCFTOOLS_NORM_CLAIRSTO.out.vcf
        .join( BCFTOOLS_NORM_CLAIRSTO.out.tbi, by: 0 )
    deepsomatic_view_ch = BCFTOOLS_NORM_DEEPSOMATIC.out.vcf
        .join( BCFTOOLS_NORM_DEEPSOMATIC.out.tbi, by: 0 )

    BCFTOOLS_VIEW_PASS_CLAIRSTO(    clairsto_view_ch,    [], [], [] )
    BCFTOOLS_VIEW_PASS_DEEPSOMATIC( deepsomatic_view_ch, [], [], [] )

    // ── bcftools isec -n+2 -c none ──────────────────────────────────────
    // BCFTOOLS_ISEC (local) takes [meta, [vcf1, vcf2], [tbi1, tbi2]] and emits
    // separate channels `.a` (first input's intersection subset) and `.b`
    // (second input's intersection subset).
    isec_input_ch = BCFTOOLS_VIEW_PASS_CLAIRSTO.out.vcf
        .join( BCFTOOLS_VIEW_PASS_CLAIRSTO.out.tbi,    by: 0 )
        .join( BCFTOOLS_VIEW_PASS_DEEPSOMATIC.out.vcf, by: 0 )
        .join( BCFTOOLS_VIEW_PASS_DEEPSOMATIC.out.tbi, by: 0 )
        .map { meta, c_vcf, c_tbi, d_vcf, d_tbi ->
            [ meta, [ c_vcf, d_vcf ], [ c_tbi, d_tbi ] ]
        }
    BCFTOOLS_ISEC( isec_input_ch )

    // ── Rename FORMAT/<source> + FORMAT/DP to caller-tagged names ───────
    // ClairS-TO uses FORMAT/AF; DeepSomatic uses FORMAT/VAF — pass per caller.
    RENAME_CLAIRSTO(    BCFTOOLS_ISEC.out.a, 'CLAIRSTO',    'AF'  )
    RENAME_DEEPSOMATIC( BCFTOOLS_ISEC.out.b, 'DEEPSOMATIC', 'VAF' )

    // ── Promote renamed FORMAT fields to INFO ───────────────────────────
    VIEW_TO_INFO_CLAIRSTO(    RENAME_CLAIRSTO.out.vcf,    ['CLAIRSTO_VAF',    'CLAIRSTO_DP']    )
    VIEW_TO_INFO_DEEPSOMATIC( RENAME_DEEPSOMATIC.out.vcf, ['DEEPSOMATIC_VAF', 'DEEPSOMATIC_DP'] )

    // ── Cross-annotate: graft DeepSomatic INFO/{VAF,DP} onto ClairS-TO ──
    // BCFTOOLS_ANNOTATE input signature (nf-core): a single tuple of 8 paths:
    //   [meta, input, index, annotations, annotations_index, columns, header_lines, rename_chrs]
    // We pass [], [], [] for the optional columns / header_lines / rename_chrs slots.
    cross_input_ch = VIEW_TO_INFO_CLAIRSTO.out.vcf
        .join( VIEW_TO_INFO_DEEPSOMATIC.out.vcf, by: 0 )
        .map { meta, c_vcf, c_tbi, d_vcf, d_tbi ->
            [ meta, c_vcf, c_tbi, d_vcf, d_tbi, [], [], [] ]
        }
    BCFTOOLS_ANNOTATE_CROSS_CALLER( cross_input_ch )

    // ── Tag INFO/CALLER=CLAIRSTO,DEEPSOMATIC ────────────────────────────
    // Local BCFTOOLS_ADD_CALLER_TAG takes [meta, vcf, tbi] + val(caller_value).
    add_caller_in_ch = BCFTOOLS_ANNOTATE_CROSS_CALLER.out.vcf
        .join( BCFTOOLS_ANNOTATE_CROSS_CALLER.out.tbi, by: 0 )
    BCFTOOLS_ADD_CALLER_TAG( add_caller_in_ch, 'CLAIRSTO,DEEPSOMATIC', [] )

    // ── Candidate-gene rescue ───────────────────────────────────────────
    // Recovers single-caller PASS variants inside rescue regions that isec -n+2
    // discarded. Runs as a PARALLEL branch, not a mix into ISEC.out.a/.b:
    // BCFTOOLS_ANNOTATE_CROSS_CALLER uses the ClairS-TO VCF as its base record
    // set and only grafts DeepSomatic INFO onto existing records, so a
    // DeepSomatic-only rescue mixed into stream b would never appear at all.
    def ch_merged_consensus = BCFTOOLS_ADD_CALLER_TAG.out.vcf

    if ( params.rescue_snv?.toString()?.toBoolean() ) {

        // Broadcastable single-item view of the rescue BED. `.first()` makes
        // this safe to reuse across every per-sample invocation below even if
        // the upstream channel is a regular (non-`Channel.value`) queue.
        def rescue_bed_path = rescue_bed.map { _m, bed -> bed }.first()

        // BCFTOOLS_VIEW (nf-core) declares 4 SEPARATE inputs — the main
        // [meta, vcf, index] tuple, then regions/targets/samples as their own
        // path channels — so the rescue BED is passed as its own positional
        // arg (regions), not folded into the main tuple.
        rescue_ct_in = BCFTOOLS_VIEW_PASS_CLAIRSTO.out.vcf
            .join( BCFTOOLS_VIEW_PASS_CLAIRSTO.out.tbi, by: 0 )
        RESCUE_VIEW_CLAIRSTO( rescue_ct_in, rescue_bed_path, [], [] )

        rescue_ds_in = BCFTOOLS_VIEW_PASS_DEEPSOMATIC.out.vcf
            .join( BCFTOOLS_VIEW_PASS_DEEPSOMATIC.out.tbi, by: 0 )
        RESCUE_VIEW_DEEPSOMATIC( rescue_ds_in, rescue_bed_path, [], [] )

        // BCFTOOLS_VIEW emits vcf/tbi as two separate named outputs (unlike
        // BCFTOOLS_ISEC, which emits [meta, vcf, tbi] as a single `.a`/`.b`
        // tuple) — rejoin them before handing off to the local rename module.
        RESCUE_RENAME_CLAIRSTO(
            RESCUE_VIEW_CLAIRSTO.out.vcf.join( RESCUE_VIEW_CLAIRSTO.out.tbi, by: 0 ),
            'CLAIRSTO', 'AF'
        )
        RESCUE_RENAME_DEEPSOMATIC(
            RESCUE_VIEW_DEEPSOMATIC.out.vcf.join( RESCUE_VIEW_DEEPSOMATIC.out.tbi, by: 0 ),
            'DEEPSOMATIC', 'VAF'
        )

        RESCUE_INFO_CLAIRSTO(    RESCUE_RENAME_CLAIRSTO.out.vcf,    ['CLAIRSTO_VAF',    'CLAIRSTO_DP']    )
        RESCUE_INFO_DEEPSOMATIC( RESCUE_RENAME_DEEPSOMATIC.out.vcf, ['DEEPSOMATIC_VAF', 'DEEPSOMATIC_DP'] )

        // Single CALLER value marks these as single-caller evidence.
        // BCFTOOLS_VIEW_TO_INFO emits vcf+tbi together as one `.vcf` tuple
        // (no separate `.tbi` output exists), so pass it straight through.
        RESCUE_TAG_CLAIRSTO(    RESCUE_INFO_CLAIRSTO.out.vcf,    'CLAIRSTO',    [] )
        RESCUE_TAG_DEEPSOMATIC( RESCUE_INFO_DEEPSOMATIC.out.vcf, 'DEEPSOMATIC', [] )

        // Transfer the matched gene symbol from the rescue BED into INFO/RESCUED.
        def rescued_hdr = file("${workflow.workDir}/rescued.hdr")
        rescued_hdr.parent.mkdirs()
        rescued_hdr.text = '##INFO=<ID=RESCUED,Number=1,Type=String,Description="Rescued single-caller variant; value is the candidate gene whose region matched">\n'

        RESCUE_FLAG_CLAIRSTO(
            RESCUE_TAG_CLAIRSTO.out.vcf
                .combine( rescue_bed_path.map { bed -> [ bed ] } )
                .map { meta, vcf, tbi, bed -> [ meta, vcf, tbi, bed, [], [], rescued_hdr, [] ] }
        )
        RESCUE_FLAG_DEEPSOMATIC(
            RESCUE_TAG_DEEPSOMATIC.out.vcf
                .combine( rescue_bed_path.map { bed -> [ bed ] } )
                .map { meta, vcf, tbi, bed -> [ meta, vcf, tbi, bed, [], [], rescued_hdr, [] ] }
        )

        // ── Normalise sample names on ALL THREE streams BEFORE concat ───────
        // `bcftools concat` (unlike `merge`) stitches records without joining
        // on sample columns, so it hard-requires every input to carry the
        // IDENTICAL sample name. ClairS-TO writes the input BAM's SM tag into
        // the sample column while DeepSomatic always writes its own fixed
        // placeholder ("Sample") — so left alone, the three concat inputs
        // (consensus, ClairS-TO rescue, DeepSomatic rescue) can carry three
        // different sample names and MERGE_SNV_CONCAT_RESCUE dies with
        // "Different sample names". Do NOT rely on the consensus already
        // being named `meta.id` — that only holds here because ClairS-TO's
        // SM tag happens to equal meta.id ("colo829") in this test dataset;
        // in production the two can differ (see the note on the final
        // BCFTOOLS_REHEADER below). The subworkflow already contains a
        // reheader-to-meta.id step, but it runs at the very end, AFTER
        // BCFTOOLS_SORT_MERGE_SNV — i.e. after this concat — so it cannot
        // repair a mismatch the concat has already failed on. Reheader here,
        // immediately before building concat_in, using the same one-line
        // samples-file pattern as the final reheader block. Do not move this
        // back downstream of the concat again.
        RESCUE_REHEADER_CONSENSUS(
            ch_merged_consensus.map { meta, vcf, tbi ->
                def sf = File.createTempFile("${meta.id}.consensus.rescue.samples", ".txt", new File("${workflow.workDir}"))
                sf.deleteOnExit()
                sf.text = "${meta.id}\n"
                [ meta, vcf, [], file(sf) ]
            },
            Channel.value( [ [:], [] ] )
        )
        ch_merged_consensus = RESCUE_REHEADER_CONSENSUS.out.vcf
            .join( RESCUE_REHEADER_CONSENSUS.out.index, by: 0 )

        RESCUE_REHEADER_CLAIRSTO(
            RESCUE_FLAG_CLAIRSTO.out.vcf
                .join( RESCUE_FLAG_CLAIRSTO.out.tbi, by: 0 )
                .map { meta, vcf, tbi ->
                    def sf = File.createTempFile("${meta.id}.clairsto.rescue.samples", ".txt", new File("${workflow.workDir}"))
                    sf.deleteOnExit()
                    sf.text = "${meta.id}\n"
                    [ meta, vcf, [], file(sf) ]
                },
            Channel.value( [ [:], [] ] )
        )

        RESCUE_REHEADER_DEEPSOMATIC(
            RESCUE_FLAG_DEEPSOMATIC.out.vcf
                .join( RESCUE_FLAG_DEEPSOMATIC.out.tbi, by: 0 )
                .map { meta, vcf, tbi ->
                    def sf = File.createTempFile("${meta.id}.deepsomatic.rescue.samples", ".txt", new File("${workflow.workDir}"))
                    sf.deleteOnExit()
                    sf.text = "${meta.id}\n"
                    [ meta, vcf, [], file(sf) ]
                },
            Channel.value( [ [:], [] ] )
        )

        // Drop rescued records whose position is ALREADY in the consensus.
        // `bcftools concat -a -D` cannot be relied on for this: it merges in
        // position order across inputs and keeps whichever duplicate it meets
        // first, which is not guaranteed to be the consensus copy. Measured on
        // real data it kept the rescue copy for 2,486 positions, downgrading
        // genuine two-caller calls to single-caller and discarding one caller's
        // VAF/DP -- which then failed FILTER_SNV's dual-VAF clause. Subtracting
        // here makes the `-D` in concat belt-and-braces rather than load-bearing.
        RESCUE_SUBTRACT_CLAIRSTO(
            RESCUE_REHEADER_CLAIRSTO.out.vcf
                .join( RESCUE_REHEADER_CLAIRSTO.out.index, by: 0 )
                .join( ch_merged_consensus, by: 0 )
                .map { meta, r_vcf, r_tbi, c_vcf, c_tbi -> [ meta, r_vcf, r_tbi, c_vcf, c_tbi ] }
        )
        RESCUE_SUBTRACT_DEEPSOMATIC(
            RESCUE_REHEADER_DEEPSOMATIC.out.vcf
                .join( RESCUE_REHEADER_DEEPSOMATIC.out.index, by: 0 )
                .join( ch_merged_consensus, by: 0 )
                .map { meta, r_vcf, r_tbi, c_vcf, c_tbi -> [ meta, r_vcf, r_tbi, c_vcf, c_tbi ] }
        )

        concat_in = ch_merged_consensus
            .join( RESCUE_SUBTRACT_CLAIRSTO.out.vcf,      by: 0 )
            .join( RESCUE_SUBTRACT_CLAIRSTO.out.tbi,      by: 0 )
            .join( RESCUE_SUBTRACT_DEEPSOMATIC.out.vcf,   by: 0 )
            .join( RESCUE_SUBTRACT_DEEPSOMATIC.out.tbi,   by: 0 )
            .map { meta, c_vcf, c_tbi, r1_vcf, r1_tbi, r2_vcf, r2_tbi ->
                [ meta, [ c_vcf, r1_vcf, r2_vcf ], [ c_tbi, r1_tbi, r2_tbi ] ]
            }
        MERGE_SNV_CONCAT_RESCUE( concat_in )
        ch_merged_consensus = MERGE_SNV_CONCAT_RESCUE.out.vcf
            .join( MERGE_SNV_CONCAT_RESCUE.out.tbi, by: 0 )
    }

    // ── Final sort (defensive; annotate usually preserves order) ────────
    sort_in_ch = ch_merged_consensus.map { meta, vcf, _tbi -> [ meta, vcf ] }
    BCFTOOLS_SORT_MERGE_SNV( sort_in_ch )

    // ── Rename sample column to meta.id ─────────────────────────────────
    // ClairS-TO / DeepSomatic emit the sample column with the caller's own
    // sample-name convention (typically the input BAM's SM tag), which can
    // differ from the sample_sheet meta.id. Left as-is, the downstream cohort
    // report groups by the VCF sample column and produces duplicate rows for
    // the same biological sample when the two names disagree. Mirrors the
    // MERGE_SV pattern (see subworkflows/local/merge_sv/main.nf around the
    // BCFTOOLS_REHEADER call): write a one-line samples file with `${meta.id}`
    // and pass it to `bcftools reheader --samples`. QC below still runs on
    // the pre-reheader VCF (sample name is orthogonal to the ground-truth
    // isec comparison); the emitted VCF uses the reheadered result.
    // NOTE: when rescue is enabled, every stream reaching this point has
    // already been reheadered to meta.id by RESCUE_REHEADER_* above (that
    // step exists because the rescue concat needs matching sample names
    // *before* it runs, not just at the end) — so this call is a harmless
    // no-op re-assertion of the same name in that path. It stays in place
    // unconditionally because it is also the ONLY normalisation step when
    // rescue is disabled (params.rescue_snv=false).
    BCFTOOLS_SORT_MERGE_SNV.out.vcf
        .map { meta, vcf ->
            def sf = File.createTempFile("${meta.id}.merged.snv.samples", ".txt", new File("${workflow.workDir}"))
            sf.deleteOnExit()
            sf.text = "${meta.id}\n"
            [ meta, vcf, [], file(sf) ]
        }
        .set { ch_reheader_snv_input }

    BCFTOOLS_REHEADER( ch_reheader_snv_input, Channel.value( [ [:], [] ] ) )

    // ── Warning-only QC: compare consensus to isec ground truth ─────────
    // QC_MERGE_SNV takes [meta, consensus_vcf, consensus_tbi, isec_0000].
    // QC uses the pre-reheader sorted VCF because the QC comparison is
    // record-level (CHROM/POS/REF/ALT) and independent of the sample column.
    qc_input_ch = BCFTOOLS_SORT_MERGE_SNV.out.vcf
        .join( BCFTOOLS_SORT_MERGE_SNV.out.tbi, by: 0 )
        .join( BCFTOOLS_ISEC.out.a.map { meta, vcf, tbi -> [ meta, vcf ] }, by: 0 )
    QC_MERGE_SNV( qc_input_ch )

    emit:
    vcf    = BCFTOOLS_REHEADER.out.vcf      // sample column renamed to meta.id
    tbi    = BCFTOOLS_REHEADER.out.index    // fresh tbi — pre-reheader tbi is
                                            // invalid once the header length
                                            // changes and BGZF offsets shift.
                                            // Written via ext.args2='--write-index=tbi'
                                            // in conf/modules.config.
    qc_log = QC_MERGE_SNV.out.log
}
