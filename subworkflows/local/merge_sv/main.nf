//
// Sort, validate, and index per-caller SV VCFs, then ensemble-merge with MINDA
//

include { BCFTOOLS_REHEADER  } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW      } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX     } from '../../../modules/nf-core/bcftools/index/main'

// options to merge the SVs
include { MINDA_MINDA as MINDA }             from '../../../modules/local/minda/minda/main'
include { SVDB_MERGE           }    from '../../../modules/nf-core/svdb/merge/main'
include { SURVIVOR_MERGE }          from '../../../modules/nf-core/survivor/merge/main'
include { JASMINESV }               from '../../../modules/nf-core/jasminesv/main'

workflow MERGE_SV {

    take:
    sv_callers_ch  // channel: [ val(meta), [ vcf, ... ] ] — one VCF list per sample
    bp_distance    // val: basepair distance tolerance for MINDA ensemble

    main:

    // Tag each VCF with its caller name by position: index 0 = severus (from samplesheet),
    // index 1 = savana (from SAVANA process). caller is embedded in meta so ext.prefix can
    // produce unique output filenames (e.g. colo829.severus.sorted.vcf.gz) and avoid
    // input/output filename collisions inside BCFTOOLS_SORT and BCFTOOLS_INDEX work dirs.
    def callers = ['severus', 'savana']
    sv_callers_ch
        .flatMap { meta, vcfs ->
            vcfs.withIndex().collect { vcf, idx ->
                [ meta + [caller: callers[idx]], vcf, [] ]
            }
        }
        .branch {
            compressed:   it[1].name.endsWith('.gz')
            uncompressed: true
        }
        .set { ch_per_vcf }

    // Only run VIEW (bgzip + validate) on files that are not already compressed
    BCFTOOLS_VIEW( ch_per_vcf.uncompressed, [], [], [] )

    // Strip the unused index placeholder from the compressed branch, then merge both paths
    ch_per_vcf.compressed
        .map { meta, vcf, _index -> [ meta, vcf ] }
        .mix( BCFTOOLS_VIEW.out.vcf )
        .branch {
            severus: it[0].caller == 'severus'
            savana:  true
        }
        .set { ch_by_caller }

    // Rename the SEVERUS sample column to meta.id.
    // bcftools reheader --samples takes a plain-text file with the new sample name(s).
    ch_by_caller.severus
        .map { meta, vcf ->
            def sf = File.createTempFile("${meta.id}.${meta.caller}.samples", ".txt", new File("${workflow.workDir}"))
            sf.deleteOnExit()
            sf.text = "${meta.id}\n"
            [ meta, vcf, [], file(sf) ]
        }
        .set { ch_reheader_input }

    BCFTOOLS_REHEADER( ch_reheader_input, Channel.value( [ [:], [] ] ) )

    // Recombine reheadered SEVERUS with SAVANA before sorting
    BCFTOOLS_REHEADER.out.vcf
        .mix( ch_by_caller.savana )
        .set { ch_sort_input }

    BCFTOOLS_SORT( ch_sort_input )
    BCFTOOLS_INDEX( BCFTOOLS_SORT.out.vcf )

    // Tap sorted+indexed channel (with caller in meta) before stripping caller.
    // Used to emit per-caller VCFs for ANNOTATIONS_SV.
    BCFTOOLS_SORT.out.vcf
        .join( BCFTOOLS_INDEX.out.tbi, by: 0 )
        .set { ch_sorted_indexed }  // [ meta+caller, vcf.gz, tbi ]

    // Strip caller from meta before groupTuple so severus and savana items —
    // which have different meta objects — are recognised as the same sample and
    // collapsed into one [ meta, [vcf1,vcf2], [tbi1,tbi2] ] tuple.
    // .map then unpacks the vcf list into the two positional paths MINDA expects.
    ch_sorted_indexed
        .map { meta, vcf, tbi ->
            [ meta.findAll { k, v -> k != 'caller' }, vcf, tbi ]
        }
        .groupTuple( by: 0 )
        .map { meta, vcfs, tbis -> [ meta, vcfs[0], vcfs[1] ] }
        .set { ch_minda }

    // 
    // MERGING  
    // 
    MINDA( ch_minda, bp_distance )

    // // svdb merge as well 
    // channel.value(params.svdb_priority) 
    //     .set{priority}
    // channel.value(params.sort_inputs.toBoolean())
    //     .set{sort_vcfs}
    // SVDB_MERGE(ch_minda, svdb_priority, sort_vcfs)
    // SURVIVOR_MERGE(ch_minda)
    // JASMINESV(ch_minda)
    
    // Per-caller emit channels (caller stripped from meta for downstream compatibility)
    ch_sorted_indexed
        .filter { meta, vcf, tbi -> meta.caller == 'severus' }
        .map    { meta, vcf, tbi -> [ meta.findAll { k, v -> k != 'caller' }, vcf, tbi ] }
        .set { ch_severus_sorted }

    ch_sorted_indexed
        .filter { meta, vcf, tbi -> meta.caller == 'savana' }
        .map    { meta, vcf, tbi -> [ meta.findAll { k, v -> k != 'caller' }, vcf, tbi ] }
        .set { ch_savana_sorted }

    emit:
    merged_sv   = MINDA.out.ensemble_vcf  // [ meta, ensemble.vcf ] — raw MINDA ensemble VCF
    severus_vcf = ch_severus_sorted       // [ meta, vcf.gz, tbi ]
    savana_vcf  = ch_savana_sorted        // [ meta, vcf.gz, tbi ]

}
