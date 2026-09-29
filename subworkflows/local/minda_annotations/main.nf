//
// Cross-reference MINDA ensemble VCF with per-caller VCFs to produce:
//   - Severus VCF annotated with CALLER + SAVANA overlap fields
//   - SAVANA-only VCF (calls absent from MINDA), sorted and indexed
//

include { MINDA_EXTRACT_ANNOTATIONS } from '../../../modules/local/minda/extract_annotations/main'
include { MINDA_FILTER_SAVANA_ONLY  } from '../../../modules/local/minda/filter_savana_only/main'
include { BCFTOOLS_ANNOTATE         } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_INDEX            } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_SORT             } from '../../../modules/nf-core/bcftools/sort/main'

workflow MINDA_ANNOTATIONS {

    take:
    minda_vcf_ch    // channel: [ val(meta), path(ensemble.vcf) ]
    severus_vcf_ch  // channel: [ val(meta), path(vcf.gz), path(tbi) ]
    savana_vcf_ch   // channel: [ val(meta), path(vcf.gz), path(tbi) ]

    main:

    // Step 1: extract SAVANA annotation data and CALLER assignments from MINDA
    // join produces: [ meta, minda_vcf, sev_vcf, sev_tbi, sav_vcf, sav_tbi ]
    minda_vcf_ch
        .join( severus_vcf_ch, by: 0 )
        .join( savana_vcf_ch,  by: 0 )
        .set { ch_extract_input }

    MINDA_EXTRACT_ANNOTATIONS( ch_extract_input )

    // Step 2: annotate the Severus VCF with CALLER + SAVANA overlap INFO fields
    // Matched by CHROM,POS,~ID in ext.args — bypasses REF/ALT allele check
    // BCFTOOLS_ANNOTATE tuple: [ meta, input, index, annotations, ann_index,
    //                            columns, header_lines, rename_chrs ]
    severus_vcf_ch
        .join( MINDA_EXTRACT_ANNOTATIONS.out.annotation_vcf, by: 0 )
        .join( MINDA_EXTRACT_ANNOTATIONS.out.header,         by: 0 )
        .map { meta, sev_vcf, sev_tbi, ann_tab, ann_tbi, header ->
            [ meta, sev_vcf, sev_tbi, ann_tab, ann_tbi, [], header, [] ]
        }
        .set { ch_annotate_input }

    BCFTOOLS_ANNOTATE( ch_annotate_input )

    // Step 3: filter SAVANA VCF to calls absent from MINDA (SAVANA-only calls)
    savana_vcf_ch
        .join( MINDA_EXTRACT_ANNOTATIONS.out.supported_ids, by: 0 )
        .join( MINDA_EXTRACT_ANNOTATIONS.out.header,        by: 0 )
        .map { meta, sav_vcf, sav_tbi, ids, header -> [ meta, sav_vcf, sav_tbi, ids, header ] }
        .set { ch_filter_input }

    MINDA_FILTER_SAVANA_ONLY( ch_filter_input )

    // Step 4: sort and index the SAVANA-only VCF
    BCFTOOLS_SORT(  MINDA_FILTER_SAVANA_ONLY.out.vcf )
    BCFTOOLS_INDEX( BCFTOOLS_SORT.out.vcf )

    emit:
    annotated_severus_vcf = BCFTOOLS_ANNOTATE.out.vcf  // [ meta, vcf.gz ] — Severus + CALLER/SAVANA INFO
    savana_only_vcf       = BCFTOOLS_SORT.out.vcf      // [ meta, vcf.gz ] — SAVANA-only calls, sorted
    savana_only_tbi       = BCFTOOLS_INDEX.out.tbi     // [ meta, tbi ]
    versions              = channel.empty()
}
