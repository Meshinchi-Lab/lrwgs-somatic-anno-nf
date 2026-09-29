//
// IGV Reports — per-sample interactive HTML for pass2 SVs (+ optionally SNVs).
//
// Architecture (per spec 2026-07-06-igv-report-expansion-design.md):
//   1. SV_TO_BEDPE converts FILTER_SV.out.vcf_pass2 (+ optionally
//      FILTER_SNV.out.vcf) to three bucket-split BEDPE files per sample.
//   2. Each bucket feeds an aliased IGVREPORTS_* invocation with bucket-
//      appropriate --flanking and --merge-overlaps=true.
//   3. All three reports show the sample BAM, user BEDs (ST17/ST19), and the
//      gnomAD SV custom track.
//
// CNAs are NOT routed through IGV — see Spec D placeholder for the future
// coverage-track-based CNA visualization workflow.
//

include { IGVREPORTS as IGVREPORTS_SMALL } from '../../../modules/local/igvreports/main'
include { IGVREPORTS as IGVREPORTS_MED   } from '../../../modules/local/igvreports/main'
include { IGVREPORTS as IGVREPORTS_LARGE } from '../../../modules/local/igvreports/main'

// bcftools query lives in its own module (bcftools 1.20 container) because
// the biocontainer has no Python; SV_TO_BEDPE runs in a Python container.
include { BCFTOOLS_QUERY_FIELDS as BCFTOOLS_QUERY_SV_FIELDS  } from '../../../modules/local/bcftools/query_fields/main'
include { BCFTOOLS_QUERY_FIELDS as BCFTOOLS_QUERY_SNV_FIELDS } from '../../../modules/local/bcftools/query_fields/main'
include { SV_TO_BEDPE                                        } from '../../../modules/local/sv_to_bedpe/main'


workflow IGV_REPORTS_SV {

    take:
    sample_bam_ch        // channel: [ meta, bam, bai ]              — haplotagged BAM per sample
    sv_vcf_pass2_ch      // channel: [ meta, vcf.gz ]                — FILTER_SV.out.vcf_pass2
    sv_tbi_pass2_ch      // channel: [ meta, vcf.gz.tbi ]            — FILTER_SV.out.tbi_pass2
    snv_vcf_ch           // channel: [ meta, vcf.gz ] or empty       — FILTER_SNV.out.vcf when enabled
    snv_tbi_ch           // channel: [ meta, vcf.gz.tbi ] or empty   — FILTER_SNV.out.tbi when enabled
    fasta_fai_ch         // channel: [ meta, fasta, fai ]
    candidate_genes_ch   // channel: value(file or [])              — canonical panel from ANNOTSV_SETUPUSERANNO.out.genes

    main:

    // ── Stage user BED files ────────────────
    // BEDs come from params.FtIncludedInSV / SVincludedInFt / AnyOverlap.
    def bed_paths = []
    if ( params.FtIncludedInSV ) bed_paths += params.FtIncludedInSV.collect { file(it, checkIfExists: true) }
    if ( params.SVincludedInFt ) bed_paths += params.SVincludedInFt.collect { file(it, checkIfExists: true) }
    if ( params.AnyOverlap     ) bed_paths += params.AnyOverlap.collect     { file(it, checkIfExists: true) }
    def bed_ch = Channel.value( bed_paths )

    // ── Stage cytoband ideogram (optional) ──
    def ch_ideogram = params.cytoband_ideo
        ? Channel.value( file(params.cytoband_ideo, checkIfExists: true) )
        : Channel.value( [] )

    // ── Stage gnomAD SV track + tbi (optional) ──
    def ch_gnomad_sv = params.gnomad_sv_track
        ? Channel.value( file(params.gnomad_sv_track, checkIfExists: true) )
        : Channel.value( [] )
    def ch_gnomad_sv_tbi = params.gnomad_sv_tbi
        ? Channel.value( file(params.gnomad_sv_tbi, checkIfExists: true) )
        : Channel.value( [] )

    // Candidate-gene panel arrives via the `candidate_genes_ch` take input —
    // canonical source is ANNOTSV_SETUPUSERANNO.out.genes (with a workflow-
    // level fallback to the raw params file when SETUPUSERANNO is skipped).
    // Passed to SV_TO_BEDPE so vcf_to_bedpe.py can prefer candidate-gene hits
    // in the BEDPE Name column when any of the record's Gene_name tokens
    // overlap the panel. Falls back to first-10 of Gene_name when empty.
    def ch_candidate_genes = candidate_genes_ch

    // ── Normalise meta on the SV pass2 channels ──
    // SV pass2 channels carry meta.variant_type='sv' (from variant_ch.svs).
    // Strip it so we can join against sample_bam_ch and snv_ch which use
    // sample-keyed meta {id: ...}.
    def strip_vt = { meta, path -> [ meta.findAll { it.key != 'variant_type' }, path ] }

    def sv_vcf_norm  = sv_vcf_pass2_ch.map( strip_vt )
    def sv_tbi_norm  = sv_tbi_pass2_ch.map( strip_vt )
    def snv_vcf_norm = snv_vcf_ch.map( strip_vt )
    def snv_tbi_norm = snv_tbi_ch.map( strip_vt )

    // ── Extract per-record fields to TSV via bcftools query ──
    // Two aliased BCFTOOLS_QUERY_FIELDS invocations run in the bcftools
    // container, each with its own -f format string. SV_TO_BEDPE then consumes
    // the TSVs in a Python-only container. The SNV path only fires when
    // snv_vcf_ch has items — when disabled (Channel.empty()), the SNV query
    // process never runs, and the join with remainder=true supplies a `[]`
    // sentinel to SV_TO_BEDPE so it skips --snv-in cleanly.
    // SV format extended with fields required to enrich the BEDPE Name column
    // — igv-reports 1.16 BedpeTable is hardcoded to 8 columns (6 coords + Name),
    // so extra annotation only surfaces if packed into Name. Order matches
    // bin/vcf_to_bedpe.py sv_row_to_bedpe() field-index expectations.
    // Fields 9+: CALLER, ACMG_class, AnnotSV_ranking_score, CytoBand, then
    // 5 name_ST17_* + 5 score_ST17_* pairs (BND/DEL/DUP/INS/INV).
    def sv_format  = '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%ALT\\t%ID\\t%INFO/Gene_name\\t%INFO/CALLER\\t%INFO/ACMG_class\\t%INFO/AnnotSV_ranking_score\\t%INFO/CytoBand\\t%INFO/name_ST17_Alterations.SV.All_BND\\t%INFO/name_ST17_Alterations.SV.All_DEL\\t%INFO/name_ST17_Alterations.SV.All_DUP\\t%INFO/name_ST17_Alterations.SV.All_INS\\t%INFO/name_ST17_Alterations.SV.All_INV\\t%INFO/score_ST17_Alterations.SV.All_BND\\t%INFO/score_ST17_Alterations.SV.All_DEL\\t%INFO/score_ST17_Alterations.SV.All_DUP\\t%INFO/score_ST17_Alterations.SV.All_INS\\t%INFO/score_ST17_Alterations.SV.All_INV\\n'
    def snv_format = '%CHROM\\t%POS\\t%ID\\t%INFO/CSQ\\n'

    sv_vcf_norm
        .join( sv_tbi_norm, by: 0 )
        .set { ch_sv_query_input }

    snv_vcf_norm
        .join( snv_tbi_norm, by: 0 )
        .set { ch_snv_query_input }

    BCFTOOLS_QUERY_SV_FIELDS ( ch_sv_query_input,  sv_format  )
    BCFTOOLS_QUERY_SNV_FIELDS( ch_snv_query_input, snv_format )

    // ── Combine SV + SNV TSVs (with SNV empty-fallback) → SV_TO_BEDPE ──
    // Single .join with remainder: true handles the SNV empty case cleanly.
    // When SNV is disabled, the right side yields null; .map substitutes [].
    BCFTOOLS_QUERY_SV_FIELDS.out.tsv
        .join( BCFTOOLS_QUERY_SNV_FIELDS.out.tsv, by: 0, remainder: true )
        .map { meta, sv_tsv, snv_tsv -> tuple( meta, sv_tsv, snv_tsv ?: [] ) }
        .set { ch_sv_to_bedpe_input }

    SV_TO_BEDPE( ch_sv_to_bedpe_input, ch_candidate_genes )

    // ── Build per-bucket IGVREPORTS input tuples ──
    // Module signature: tuple val(meta), path(sites), path(tracks), path(tracks_indices)
    //   tracks         = [ BAM, *BED files ]
    //   tracks_indices = [ BAI ]  — BEDPE has no companion index; TBI list stays empty
    def sample_bam_norm = sample_bam_ch
        .map { meta, bam, bai -> [ meta.findAll { it.key != 'variant_type' }, bam, bai ] }

    sample_bam_norm
        .join( SV_TO_BEDPE.out.small, by: 0 )
        .combine( bed_ch )
        .map { it ->
            def meta = it[0]
            def bam  = it[1]
            def bai  = it[2]
            def sites = it[3]
            def beds  = it.size() > 4 ? it[4..-1] : []
            tuple( meta, sites, [ bam ] + beds, [ bai ] )
        }
        .set { ch_igv_small }

    sample_bam_norm
        .join( SV_TO_BEDPE.out.medium, by: 0 )
        .combine( bed_ch )
        .map { it ->
            def meta = it[0]
            def bam  = it[1]
            def bai  = it[2]
            def sites = it[3]
            def beds  = it.size() > 4 ? it[4..-1] : []
            tuple( meta, sites, [ bam ] + beds, [ bai ] )
        }
        .set { ch_igv_med }

    sample_bam_norm
        .join( SV_TO_BEDPE.out.large, by: 0 )
        .combine( bed_ch )
        .map { it ->
            def meta = it[0]
            def bam  = it[1]
            def bai  = it[2]
            def sites = it[3]
            def beds  = it.size() > 4 ? it[4..-1] : []
            tuple( meta, sites, [ bam ] + beds, [ bai ] )
        }
        .set { ch_igv_large }

    IGVREPORTS_SMALL( ch_igv_small, fasta_fai_ch, ch_ideogram, ch_gnomad_sv, ch_gnomad_sv_tbi )
    IGVREPORTS_MED  ( ch_igv_med,   fasta_fai_ch, ch_ideogram, ch_gnomad_sv, ch_gnomad_sv_tbi )
    IGVREPORTS_LARGE( ch_igv_large, fasta_fai_ch, ch_ideogram, ch_gnomad_sv, ch_gnomad_sv_tbi )

    emit:
    report        = IGVREPORTS_SMALL.out.report
                        .mix( IGVREPORTS_MED.out.report )
                        .mix( IGVREPORTS_LARGE.out.report )
    report_small  = IGVREPORTS_SMALL.out.report
    report_medium = IGVREPORTS_MED.out.report
    report_large  = IGVREPORTS_LARGE.out.report
    versions      = IGVREPORTS_SMALL.out.versions
                        .mix( IGVREPORTS_MED.out.versions )
                        .mix( IGVREPORTS_LARGE.out.versions )
}
