//
// Deduplicate SAVANA same-position INS records before MERGE_SV/MINDA sees them.
//
// SAVANA emits multiple INS calls at the same breakpoint with different inserted-
// sequence lengths. AnnotSV's primary key (chrom_start_end_svtype_n) collapses
// these to one AnnotSV_ID, and variantconvert then aborts with "Each variant is
// assumed to only have one single line of 'full' annotation" (variantconvert
// issue #20, closed by maintainer as by-design).
//
// Strategy: split the SAVANA VCF by SVTYPE so the dedup script only walks INS
// records (typically a small fraction of the total). Non-INS records pass
// through untouched. Recombine + sort + index.
//
// We use bcftools view (not vembrane filter) for the SVTYPE split: vembrane
// special-cases BND records by keeping mate pairs adjacent in the output, which
// breaks positional sort when the mates live at distant genomic coordinates.
// bcftools view preserves input order strictly and emits bgzipped output
// directly, so no separate bgzip step is needed.
//
// Selection cascade for which INS to keep per (CHROM, POS):
//   highest QUAL  →  highest TUMOUR_READ_SUPPORT  →  longest SVLEN
//

include { BCFTOOLS_SORT   as SAVANA_DEDUP_SORT_INPUT } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW   as BCFTOOLS_VIEW_INS       } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW   as BCFTOOLS_VIEW_NON_INS   } from '../../../modules/nf-core/bcftools/view/main'
include { DEDUPLICATE_INS                            } from '../../../modules/local/deduplicate_ins/main'
include { BCFTOOLS_CONCAT as SAVANA_DEDUP_CONCAT     } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT   as SAVANA_DEDUP_SORT       } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX  as SAVANA_DEDUP_INDEX      } from '../../../modules/nf-core/bcftools/index/main'

workflow DEDUPLICATE_SAVANA {

    take:
    savana_vcf_ch  // channel: [ val(meta), path(vcf) ] — raw SAVANA SV VCF

    main:

    // SAVANA can emit BND mate pairs in mate-pair order rather than strict
    // positional order, so the raw VCF is sometimes locally unsorted. Pre-sort
    // once at the top so both filter branches inherit sorted input — required
    // for tabix-indexing later in BCFTOOLS_CONCAT.
    SAVANA_DEDUP_SORT_INPUT( savana_vcf_ch )

    // Wrap the sorted VCF with an empty index slot for BCFTOOLS_VIEW's input contract.
    SAVANA_DEDUP_SORT_INPUT.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_view_input }

    // Split SVTYPE == INS vs everything else using bcftools view. The SVTYPE
    // expressions ride in ext.args (modules.config); only bgzipped output here.
    BCFTOOLS_VIEW_INS    ( ch_view_input, [], [], [] )
    BCFTOOLS_VIEW_NON_INS( ch_view_input, [], [], [] )

    DEDUPLICATE_INS( BCFTOOLS_VIEW_INS.out.vcf )

    // Recombine: deduped INS + bgzipped non-INS.
    DEDUPLICATE_INS.out.vcf
        .join( BCFTOOLS_VIEW_NON_INS.out.vcf, by: 0 )
        .map { meta, ins_vcf, non_ins_vcf -> [ meta, [ ins_vcf, non_ins_vcf ], [] ] }
        .set { ch_concat_input }

    SAVANA_DEDUP_CONCAT( ch_concat_input )
    SAVANA_DEDUP_SORT  ( SAVANA_DEDUP_CONCAT.out.vcf )
    SAVANA_DEDUP_INDEX ( SAVANA_DEDUP_SORT.out.vcf )

    emit:
    vcf      = SAVANA_DEDUP_SORT.out.vcf    // [ meta, *.savana.deduped.sorted.vcf.gz ]
    tbi      = SAVANA_DEDUP_INDEX.out.tbi   // [ meta, *.savana.deduped.sorted.vcf.gz.tbi ]
    versions = channel.empty()
}
