//
// COHORT_MERGE_SV -- cluster structural variants ACROSS cohort samples to identify recurrent variants.
//
// Runs AFTER SV_REPORT_INDEX: driver_tier is computed in the Quarto layer by assign_driver_tier() 
// The BED built from the tiered TSV defines REGIONS to investigate, but the records searched are the PRE-FILTER, VEP-annotated VCFs. 
// to find sub-threshold evidence in samples where the same SV was called but filtered out. 
// Region-restriction collapses this to the loci  where at least one sample already has a high-confidence call.
//
//

include { COLLECT_COORDINATES                          } from '../../../modules/local/collect_coordinates/main'
include { BCFTOOLS_ANNOTATE     as COHORT_SV_SET_ID    } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW         as COHORT_SV_RESTRICT  } from '../../../modules/nf-core/bcftools/view/main'
include { JASMINESV             as JASMINE_COHORT      } from '../../../modules/nf-core/jasminesv/main'
include { BCFTOOLS_SORT         as COHORT_SV_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_QUERY_FIELDS as COHORT_SV_RECORDS   } from '../../../modules/local/bcftools/query_fields/main'
include { JASMINE_CLUSTER_MAP                          } from '../../../modules/local/jasmine_cluster_map/main'

workflow COHORT_MERGE_SV {
    take:
    sv_tsv        // channel: [ val(meta), path(cohort_sv_pass2.tsv) ]  from SV_REPORT_INDEX
    sv_vcfs       // channel: [ val(meta), path(vcf), path(tbi) ]  ANNOTATIONS_SV.out.vep_vcf
                  //          joined to vep_tbi -- PRE-FILTER, VEP/gnomAD/ClinVar annotated
    fasta         // channel: [ val(meta), path(fasta) ]
    fasta_fai     // channel: [ val(meta), path(fai) ]

    main:

    // Fail fast: a pad smaller than the clustering distance silently excludes the very
    // records the merge exists to find, rather than producing an error.
    if ( params.cohort_sv_region_pad < params.cohort_sv_max_dist ) {
        error "cohort_sv_region_pad (${params.cohort_sv_region_pad}) must be >= " +
              "cohort_sv_max_dist (${params.cohort_sv_max_dist}); a smaller pad " +
              "excludes offset breakpoints before Jasmine can cluster them."
    }

    // Fail fast: ANNOTATIONS_SV emits vep_vcf as channel.empty() when VEP is skipped.
    // Without this the subworkflow would silently produce nothing, which looks identical
    // to "no recurrent SVs found".
    // NB: sv_anno_vep is the STRING 'true' in nextflow.config, so a bare truth test is
    // wrong -- 'false' is also a truthy String. Mirrors annotations_sv/main.nf.
    def run_vep = params.sv_anno_vep ? params.sv_anno_vep.toString().toBoolean() : false
    if ( !run_vep ) {
        error "COHORT_MERGE_SV requires the VEP-annotated pre-filter VCF " +
              "(ANNOTATIONS_SV.out.vep_vcf), but sv_anno_vep is false so that channel " +
              "is empty. Set --sv_anno_vep true, or --run_cohort_merge_sv false."
    }

    // -- Regions: tiered coordinates from the report TSV, padded --------------
    COLLECT_COORDINATES( sv_tsv, params.cohort_sv_region_pad, '1,2,3' )

    // -- Uniquify IDs per sample so Jasmine's IDLIST is unambiguous ------------
    // Replaces the reference repo's make_unique_ids.py with a native bcftools call.
    // Two samples can both emit `severus_BND_ST2671_1`; prefixing with meta.id makes the
    // cluster map joinable back to (SAMPLE, ID) in the report TSV.
    COHORT_SV_SET_ID(
        sv_vcfs.map { meta, vcf, tbi -> [ meta, vcf, tbi, [], [], [], [], [] ] }
    )

    // -- Restrict to the padded tiered regions --------------------------------
    COHORT_SV_RESTRICT(
        COHORT_SV_SET_ID.out.vcf
            .join( COHORT_SV_SET_ID.out.tbi, by: 0 ),
        COLLECT_COORDINATES.out.bed.map { _m, bed -> bed },
        [],
        []
    )

    // -- Cluster across samples ------------------------------------------------
    // Jasmine takes ALL sample VCFs as ONE input tuple, so collect() first. The module
    // declares path(vcfs, arity:'1..*') plus bams/bais/sample_dists, unused here.
    ch_collected_vcfs = COHORT_SV_RESTRICT.out.vcf
        .map { _meta, vcf -> vcf }
        .collect()

    // cohort merging is undefined for a single sample.
    // Filtering the collected channel empties it, so JASMINE_COHORT, COHORT_SV_SORT and
    // JASMINE_CLUSTER_MAP never execute if the inputs never arrive.  also covers the case where the samplesheet has
    // N>=2 but filtering leaves only one sample with records.
    ch_jasmine_in = ch_collected_vcfs
        .filter { vcfs ->
            def n = vcfs instanceof List ? vcfs.size() : 1
            if ( n < 2 ) {
                log.warn "COHORT_MERGE_SV: only ${n} sample(s) have records in the tiered " +
                         "regions; cohort SV merging requires >= 2. Skipping Jasmine and all " +
                         "downstream clustering. The regions BED and per-sample restricted " +
                         "VCFs are still published."
                return false
            }
            return true
        }
        .map { vcfs -> [ [ id: 'cohort' ], vcfs, [], [], [] ] }

    JASMINE_COHORT( ch_jasmine_in, fasta, fasta_fai, [] )

    COHORT_SV_SORT( JASMINE_COHORT.out.vcf )

    // -- Per-record annotations for RESCUED rows -------------------------------
    // Jasmine emits ONE representative record per cluster, so member-level annotations are
    // not recoverable from the merged VCF. Rescued rows therefore take theirs from the
    // region-restricted per-sample VCFs, queried natively. collectFile concatenates the
    // per-sample projections into one cohort-wide table -- no custom script.
    ch_records = COHORT_SV_RECORDS(
            COHORT_SV_RESTRICT.out.vcf.join( COHORT_SV_RESTRICT.out.tbi, by: 0 ),
            '[%SAMPLE]\\t%ID\\t%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN' +
            '\\t%INFO/Gene_name\\t%INFO/gnomAD_SV_AF\\t%INFO/ClinVar_SV_CLNSIG\\t%INFO/CALLER\\n'
        )
        .tsv
        .map { _meta, tsv -> tsv }
        .collectFile( name: 'cohort_sv_records.tsv', sort: true )

    // -- Cluster map -> outer-joined recurrence TSV ----------------------------
    // Three-way join: cluster map (cluster<->member) x record table (member->annotations)
    // x tiered TSV (member->driver_tier). OUTER on the tiered side so sub-threshold
    // cluster members survive as RESCUED rows.
    JASMINE_CLUSTER_MAP(
        COHORT_SV_SORT.out.vcf.join( COHORT_SV_SORT.out.tbi, by: 0 ),
        sv_tsv.map { _m, tsv -> tsv },
        ch_records,
        params.cohort_sv_min_recurrence,
        params.cohort_sv_max_dist
    )

    emit:
    bed         = COLLECT_COORDINATES.out.bed
    merged_vcf  = COHORT_SV_SORT.out.vcf
    recurrent   = JASMINE_CLUSTER_MAP.out.tsv
    cluster_map = JASMINE_CLUSTER_MAP.out.map
    versions    = channel.empty()
}
