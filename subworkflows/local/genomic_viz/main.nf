// subworkflows/local/genomic_viz/main.nf
//
// GENOMIC_VIZ -- static cohort circos (+ one per sample) and a cohort oncoprint from the
// report's pass-1 cohort TSVs, using Circos and maftools (standard tools, standard inputs).
//
// Runs between the two SV_REPORT_INDEX passes, in the slot COHORT_MERGE_SV occupies, because
// driver_tier and the recurrence map exist only after pass 1. Pass 2 stages everything this
// emits under viz/ and renders the "Cohort figures" page.
//
include { COHORT_TO_CIRCOS  } from '../../../modules/local/cohort_to_circos/main'
include { COHORT_TO_MAF     } from '../../../modules/local/cohort_to_maf/main'
include { CIRCOS            } from '../../../modules/local/circos/main'
include { MAFTOOLS_ONCOPLOT } from '../../../modules/local/maftools/oncoplot/main'

workflow GENOMIC_VIZ {
    take:
    sv_tsv           // channel: [ val(meta), path(cohort_sv_pass2.tsv) ]; meta must be [ id: 'cohort' ]:
                     //   meta.id names the oncoprint files (cohort.oncoprint.*) the report reads
    cna_tsv          // value channel: path(cohort_cna_pass2.tsv) or []
    snv_tsv          // value channel: path(cohort_snv.tsv) or []
    recurrence_map   // value channel: path(cohort_sv_recurrent.map.tsv) or []
    candidate_genes  // value channel: path(candidate_genes.txt) or []
    circos_conf      // value channel: list of the assets/genomic_viz/circos/*.conf files (files, not the dir, so -resume sees edits)
    vep_lookup       // value channel: path(assets/genomic_viz/vep_consequence_to_maf.tsv)
    top_genes        // value channel: int

    main:

    // -- Circos: convert, then fan one CIRCOS task out per track set --------------
    // The track sets are 'cohort' plus one per sample. 'cohort' is reserved: COHORT_TO_CIRCOS
    // rejects a sample of that name, so the cohort set never collides with a sample set.
    COHORT_TO_CIRCOS( sv_tsv, cna_tsv, snv_tsv, recurrence_map )

    // Each converter output is one [ meta, files ] tuple (a bare path when the cohort is the
    // only set): fan it out and key each file by set name. The suffix is cut from the END of
    // the name (String.minus would remove its first occurrence anywhere), and the joins fail
    // on a mismatch: a set missing one of its three files is a converter bug, not a figure to
    // drop silently. Local closures are invoked with .call(): the strict-syntax linter
    // (Nextflow 25.10) reads a bare name(args) as a call to an undefined named function.
    def by_set = { ch, suffix ->
        ch.transpose()
          .map { _meta, f -> [ f.name.substring(0, f.name.length() - suffix.length()), f ] }
    }
    ch_sets = by_set.call( COHORT_TO_CIRCOS.out.links,  '.links.txt'  )
        .join( by_set.call( COHORT_TO_CIRCOS.out.cna,    '.cna.txt'    ), by: 0, failOnMismatch: true )
        .join( by_set.call( COHORT_TO_CIRCOS.out.labels, '.labels.txt' ), by: 0, failOnMismatch: true )
        .map { set, links, cna, labels ->
            [ [ id: set, layers: set == 'cohort' ? 6 : 1 ], links, cna, labels ]
        }

    CIRCOS( ch_sets, circos_conf )

    // -- Oncoprint ---------------------------------------------------------------
    COHORT_TO_MAF( sv_tsv, cna_tsv, snv_tsv, candidate_genes, vep_lookup, top_genes )

    MAFTOOLS_ONCOPLOT(
        COHORT_TO_MAF.out.maf
            .join( COHORT_TO_MAF.out.cn_table, by: 0 )
            .join( COHORT_TO_MAF.out.genes,    by: 0 )
            .join( COHORT_TO_MAF.out.clinical, by: 0 )
    )

    // -- Everything pass 2 stages under viz/ -------------------------------------
    // Every source is [ meta, path-or-paths ]: mix the tuples, then drop meta once. The
    // parentheses around the ternary are required by the 25.10 parser.
    ch_viz_files = CIRCOS.out.png
        .mix( CIRCOS.out.svg,
              COHORT_TO_CIRCOS.out.links, COHORT_TO_CIRCOS.out.cna, COHORT_TO_CIRCOS.out.labels, COHORT_TO_CIRCOS.out.counts,
              COHORT_TO_MAF.out.maf, COHORT_TO_MAF.out.cn_table, COHORT_TO_MAF.out.clinical,
              COHORT_TO_MAF.out.genes, COHORT_TO_MAF.out.unmapped, COHORT_TO_MAF.out.counts,
              MAFTOOLS_ONCOPLOT.out.png, MAFTOOLS_ONCOPLOT.out.svg, MAFTOOLS_ONCOPLOT.out.matrix,
              MAFTOOLS_ONCOPLOT.out.summaries, MAFTOOLS_ONCOPLOT.out.no_data )
        .flatMap { _meta, files -> (files instanceof List ? files : [ files ]) }
        .collect()

    emit:
    viz_files     = ch_viz_files                 // [ png, svg, txt, ... ] -- one list; emits nothing when
                                                 //   sv_tsv is empty, so callers apply .ifEmpty([])
    circos_png    = CIRCOS.out.png               // [ meta, png ] per set
    circos_svg    = CIRCOS.out.svg
    oncoprint_png = MAFTOOLS_ONCOPLOT.out.png
    oncoprint_svg = MAFTOOLS_ONCOPLOT.out.svg
    counts        = COHORT_TO_CIRCOS.out.counts
    versions      = channel.empty()
}
