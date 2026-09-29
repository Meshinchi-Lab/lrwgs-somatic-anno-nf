// subworkflows/local/annotations_snv/main.nf
//
// Annotate the SNV consensus VCF with ClinVar (top-level INFO via bcftools annotate),
// VEP (CSQ via the nf-core ensemblvep module), and CIViC (CIVIC INFO via bcftools annotate).
//
// Order: ClinVar → VEP → CIViC → sort.
//   - ClinVar runs FIRST: streaming and cheap; VEP is heavyweight.
//   - CIViC runs LAST, after VEP, so VEP can never clobber the CIVIC field.
//     SETUPCIVIC renames CIViC's own CSQ to CIVIC, so it cannot collide by name
//     with the CSQ that VEP writes — running CIViC last just removes the
//     question entirely, and bcftools annotate never strips existing INFO.
//     Safe because BCFTOOLS_NORM upstream in MERGE_SNV has already split
//     multi-allelics, so VEP's --pick_allele_gene can't strand a record that
//     CIViC still needs to annotate.
//
// The CIViC reference is annotated with `bcftools annotate` against the
// GRCh38-lifted VCF produced by PREPARE_REFERENCES. CIViC publishes GRCh37
// coordinates only, so annotating a GRCh38 callset against the raw release
// matched nothing and CIVIC came back empty for every record.

include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_CLINVAR_SNV } from '../../../modules/nf-core/bcftools/annotate/main'
include { ENSEMBLVEP_VEP    as ENSEMBLVEP_VEP_SNV            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_CIVIC_SNV   } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_SORT     as BCFTOOLS_SORT_ANNOTATIONS_SNV } from '../../../modules/nf-core/bcftools/sort/main'

workflow ANNOTATIONS_SNV {
    take:
    snv_vcf        // channel: [ val(meta), path(consensus.snv.vcf.gz) ]
    snv_tbi        // channel: [ val(meta), path(consensus.snv.vcf.gz.tbi) ]
    clinvar_vcf    // value:   path(clinvar.norm.vcf.gz)
    clinvar_tbi    // value:   path(clinvar.norm.vcf.gz.tbi)
    civic_vcf      // value:   path(civic.GRCh38.vcf.gz) — lifted by PREPARE_REFERENCES
    civic_tbi      // value:   path(civic.GRCh38.vcf.gz.tbi)
    vep_cache      // value:   [ [:], path(vep_cache_dir) or [] ]
    fasta          // value:   [ [meta], path(fasta) ]

    main:

    //  ClinVar annotation via bcftools annotate 
    // BCFTOOLS_ANNOTATE input:
    //   [meta, input, index, annotations, annotations_index, columns, header_lines, rename_chrs]
    // The -c column-spec lives in ext.args (modules.config); columns_file/header_lines/rename_chrs are [].
    clinvar_in_ch = snv_vcf
        .join( snv_tbi, by: 0 )
        .combine( clinvar_vcf.map { v -> [ v ] } )
        .combine( clinvar_tbi.map { t -> [ t ] } )
        .map { meta, vcf, tbi, cv, ct -> [ meta, vcf, tbi, cv, ct, [], [], [] ] }
    BCFTOOLS_ANNOTATE_CLINVAR_SNV( clinvar_in_ch )

    //  VEP consequence + IMPACT + gnomAD AF in CSQ 
    // ENSEMBLVEP_VEP first input: [ meta, vcf, custom_extra_files ].
    vep_in_ch = BCFTOOLS_ANNOTATE_CLINVAR_SNV.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
    ENSEMBLVEP_VEP_SNV(
        vep_in_ch,
        params.genome,
        params.species,
        params.cache_version,
        vep_cache,
        fasta,
        []
    )

    //  CIViC somatic mutation database
    // bcftools annotate against the GRCh38-lifted CIViC reference from PREPARE_REFERENCES. 
    civic_in_ch = ENSEMBLVEP_VEP_SNV.out.vcf
        .join( ENSEMBLVEP_VEP_SNV.out.tbi, by: 0 )
        .combine( civic_vcf.map { v -> [ v ] } )
        .combine( civic_tbi.map { t -> [ t ] } )
        .map { meta, vcf, tbi, cv, ct -> [ meta, vcf, tbi, cv, ct, [], [], [] ] }
    BCFTOOLS_ANNOTATE_CIVIC_SNV( civic_in_ch )

    //  Sort + index final annotated VCF 
    BCFTOOLS_SORT_ANNOTATIONS_SNV( BCFTOOLS_ANNOTATE_CIVIC_SNV.out.vcf )

    emit:
    vcf = BCFTOOLS_SORT_ANNOTATIONS_SNV.out.vcf
    tbi = BCFTOOLS_SORT_ANNOTATIONS_SNV.out.tbi
}
