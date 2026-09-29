// subworkflows/local/prepare_references/main.nf
//
// Stage every annotation reference the pipeline needs, once per run, so the
// variant-flow subworkflows receive ready-to-use channels instead of each
// re-deriving reference state.
//
// Covers:
//   VEP cache            — use params.vep_cache, else ENSEMBLVEP_DOWNLOAD, else skip
//   AnnotSV annotations  — use params.annotsv_annotations, else ANNOTSV_INSTALLANNOTATIONS
//   AnnotSV user BEDs    — ANNOTSV_SETUPUSERANNO (ST17 SV + ST19 CNA in one call)
//   ClinVar              — SETUPCLINVAR (chr-rename, norm, sort, index)
//   CIViC                — SETUPCIVIC → sanitize → GRCh37→GRCh38 liftover → index
//
// The CIViC chain exists because CIViC publishes GRCh37 coordinates only
// (`civicpy create-vcf` has no build flag). Annotating a GRCh38 callset against
// it silently matched nothing: BRAF V600E sits at chr7:140453136 in CIViC but
// chr7:140753336 in GRCh38, so bcftools annotate returned an empty CIVIC field
// for every record. picard LiftoverVcf with the hg19→hg38 chain fixes the
// coordinates; BCFTOOLS_SANITIZE_VCF first makes the file acceptable to
// htsjdk, which is far stricter than htslib about the VCF spec.

include { ENSEMBLVEP_DOWNLOAD                    } from '../../../modules/nf-core/ensemblvep/download/main'
include { ANNOTSV_INSTALLANNOTATIONS             } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_SETUPUSERANNO                  } from '../../../modules/local/annotsv/setupuserannotations/main'
include { BUILD_RESCUE_BED                       } from '../../../modules/local/build_rescue_bed/main'
include { SETUPCLINVAR                           } from '../../../modules/local/setupclinvar/main'
include { SETUPCLINVAR as SETUPCIVIC             } from '../../../modules/local/setupclinvar/main'
include { BCFTOOLS_SANITIZE_VCF                  } from '../../../modules/local/bcftools/sanitize_vcf/main'
include { PICARD_CREATESEQUENCEDICTIONARY        } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_LIFTOVERVCF as LIFTOVER_CIVIC   } from '../../../modules/nf-core/picard/liftovervcf/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CIVIC } from '../../../modules/nf-core/bcftools/index/main'

workflow PREPARE_REFERENCES {
    take:
    fasta          // channel: [ val(meta), path(fasta) ] — liftover target reference
    prepare_snv    // val: stage the ClinVar + CIViC SNV references (params.run_snv_anno)

    main:

    // ── VEP cache ───────────────────────────────────────────────────────
    //   local cache provided                                  → use it
    //   no cache, sv_anno_vep && include_cache                → download (published for reuse)
    //   otherwise                                             → empty tuple; VEP runs --offline or is skipped
    def ch_vep_cache
    if ( params.vep_cache ) {
        ch_vep_cache = Channel.value( [ [:], file(params.vep_cache, checkIfExists: true) ] )
    } else if ( params.sv_anno_vep?.toString()?.toBoolean() && params.include_cache?.toString()?.toBoolean() ) {
        ENSEMBLVEP_DOWNLOAD(
            Channel.value( [ [:], params.genome, params.species, params.cache_version ] ),
            Channel.value( true )
        )
        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
    } else {
        ch_vep_cache = Channel.value( [ [:], [] ] )
    }

    // ── AnnotSV annotation bundle ───────────────────────────────────────
    def ch_annotsv_base
    if ( params.annotsv_annotations ) {
        ch_annotsv_base = Channel.value([ [:], file(params.annotsv_annotations, checkIfExists: true) ])
    } else {
        ANNOTSV_INSTALLANNOTATIONS()
        ch_annotsv_base = ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { ann -> [ [:], ann ] }
            .first()
    }

    // ── AnnotSV user annotations (ST17 SV + ST19 CNA BEDs in one call) ──
    def ft_beds_combined  = (params.FtIncludedInSV ?: []) + (params.cna_FtIncludedInSV ?: [])
    def sv_beds_combined  = (params.SVincludedInFt ?: []) + (params.cna_SVincludedInFt ?: [])
    def any_beds_combined = (params.AnyOverlap     ?: []) + (params.cna_AnyOverlap     ?: [])
    def has_user_beds     = ft_beds_combined || sv_beds_combined || any_beds_combined

    // Canonical candidate-gene panel. When SETUPUSERANNO runs it re-emits this
    // file through `.out.genes`, which becomes the single downstream source;
    // when it doesn't run, fall back to the same raw file so the channel shape
    // is identical either way.
    def cand_genes_file = params.candidateGenesFile
        ? file(params.candidateGenesFile, checkIfExists: true)
        : []

    def ch_annotsv_ann
    def ch_candidate_genes_file
    if ( has_user_beds ) {
        def ft_files  = ft_beds_combined.collect  { file(it, checkIfExists: true) }
        def sv_files  = sv_beds_combined.collect  { file(it, checkIfExists: true) }
        def any_files = any_beds_combined.collect { file(it, checkIfExists: true) }

        // Warn-only header check. Without the column-header line AnnotSV emits
        // INFO keys as '_<stem>' / '_<stem>.1' instead of 'name_<stem>' /
        // 'score_<stem>'; ANNOTSV_ANNOTSV reheaders them downstream, but
        // KnotAnnotSV tooltips lose the BED attribution at the source.
        (ft_files + sv_files + any_files).each { bed ->
            def first_real = bed.readLines().find { it.trim() }
            def looks_like_header = first_real ==~ /(?i)^#\s*chrom\b.*/
            if ( !looks_like_header ) {
                log.warn """[PREPARE_REFERENCES] BED file '${bed}' is missing the column-header line.
  AnnotSV will emit INFO keys as '_<stem>' / '_<stem>.1' / '_<stem>.2' instead of
  'name_<stem>' / 'score_<stem>' / 'strand_<stem>'.
  To fix at the source, prepend this line to the BED:
      #chrom\tchromStart\tchromEnd\tname\tscore\tstrand
""".stripIndent()
            }
        }

        ANNOTSV_SETUPUSERANNO(
            ch_annotsv_base,
            Channel.value( ft_files ),
            Channel.value( sv_files ),
            Channel.value( any_files ),
            params.genome,
            cand_genes_file
        )
        ch_annotsv_ann          = ANNOTSV_SETUPUSERANNO.out.annotations
        ch_candidate_genes_file = ANNOTSV_SETUPUSERANNO.out.genes.ifEmpty( [] )
    } else {
        // No user BEDs — skip SETUPUSERANNO rather than re-copy the whole
        // annotations dir just to stage a small text file.
        ch_annotsv_ann          = ch_annotsv_base
        ch_candidate_genes_file = Channel.value( cand_genes_file )
    }

    // Tuple form required by ANNOTATIONS_SV / ANNOTATIONS_CNA `take:` inputs.
    // Groovy treats `[]` as falsy, so the ternary handles the absent panel.
    def ch_candidate_genes = ch_candidate_genes_file
        .map { f -> f ? [ [ id: f.simpleName ], f ] : [ [:], [] ] }

    // ── SNV references: ClinVar + CIViC ─────────────────────────────────
    def ch_clinvar_vcf = Channel.value( [] )
    def ch_clinvar_tbi = Channel.value( [] )
    def ch_civic_vcf   = Channel.value( [] )
    def ch_civic_tbi   = Channel.value( [] )

    if ( prepare_snv ) {

        SETUPCLINVAR(
            file(params.clinvar_snv_vcf, checkIfExists: true),
            file("${params.clinvar_snv_vcf}.tbi", checkIfExists: true),
            []
        )
        ch_clinvar_vcf = SETUPCLINVAR.out.vcf
        ch_clinvar_tbi = SETUPCLINVAR.out.tbi

        // CIViC step 1 — rename INFO/CSQ to INFO/CIVIC so it cannot collide
        // with the CSQ that VEP writes, and chr-prefix the contigs (CIViC ships
        // NCBI-style '7'; both the chain file and the callset use 'chr7').
        // SETUPCLINVAR's third input is a plain `path`, so materialise the
        // two-column rename map once here rather than routing it via a channel.
        def civic_rename = file("${workflow.workDir}/civic.rename_annots.txt")
        civic_rename.parent.mkdirs()
        civic_rename.text = "INFO/CSQ INFO/CIVIC\n"

        SETUPCIVIC(
            file(params.civic_snv_vcf, checkIfExists: true),
            file("${params.civic_snv_vcf}.tbi", checkIfExists: true),
            civic_rename
        )

        // CIViC step 2 — make it acceptable to htsjdk (see module header).
        // SETUPCLINVAR emits bare `path` channels, so attach a meta here; both
        // are single-element channels, hence combine rather than join.
        SETUPCIVIC.out.vcf
            .combine( SETUPCIVIC.out.tbi )
            .map { vcf, tbi -> [ [ id: 'civic' ], vcf, tbi ] }
            .set { ch_civic_staged }

        BCFTOOLS_SANITIZE_VCF( ch_civic_staged )

        // CIViC step 3 — GRCh37 → GRCh38. LiftoverVcf needs a sequence
        // dictionary whose basename matches the FASTA (ext.prefix in
        // modules.config enforces that).
        PICARD_CREATESEQUENCEDICTIONARY( fasta )

        LIFTOVER_CIVIC(
            BCFTOOLS_SANITIZE_VCF.out.vcf,
            PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
            fasta,
            Channel.value( [ [:], file(params.civic_chain, checkIfExists: true) ] )
        )

        // LiftoverVcf writes coordinate-sorted output but no index.
        BCFTOOLS_INDEX_CIVIC( LIFTOVER_CIVIC.out.vcf_lifted )

        // Drop the meta so CIViC matches the bare-path shape of clinvar_vcf /
        // clinvar_tbi that ANNOTATIONS_SNV already consumes.
        ch_civic_vcf = LIFTOVER_CIVIC.out.vcf_lifted.map { _meta, vcf -> vcf }
        ch_civic_tbi = BCFTOOLS_INDEX_CIVIC.out.tbi.map  { _meta, tbi -> tbi }
    }

    // ── SNV rescue regions ──────────────────────────────────────────────
    def ch_rescue_bed = Channel.value( [ [:], [] ] )
    if ( params.rescue_snv?.toString()?.toBoolean() ) {
        def rescue_genes = params.rescue_genes_file
            ? file(params.rescue_genes_file, checkIfExists: true)
            : ( params.candidateGenesFile ? file(params.candidateGenesFile, checkIfExists: true) : [] )
        if ( rescue_genes ) {
            BUILD_RESCUE_BED( ch_annotsv_ann, rescue_genes, params.genome, params.rescue_promoter_pad )
            ch_rescue_bed = BUILD_RESCUE_BED.out.bed
        } else {
            log.warn "[PREPARE_REFERENCES] rescue_snv is true but no gene list is set; rescue disabled."
        }
    }

    emit:
    vep_cache            = ch_vep_cache              // [ meta, cache_dir ] or [ [:], [] ]
    annotsv_annotations  = ch_annotsv_ann            // [ meta, annotations_dir ]
    candidate_genes      = ch_candidate_genes        // [ meta, panel ] or [ [:], [] ]
    candidate_genes_file = ch_candidate_genes_file   // raw path or []
    clinvar_vcf          = ch_clinvar_vcf            // [ meta, vcf ]
    clinvar_tbi          = ch_clinvar_tbi            // [ meta, tbi ]
    civic_vcf            = ch_civic_vcf              // [ meta, vcf ] — GRCh38, lifted
    civic_tbi            = ch_civic_tbi              // [ meta, tbi ]
    rescue_bed           = ch_rescue_bed              // [ meta, bed ] or [ [:], [] ]
}
