process IGVREPORTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7e72e81_0':
        'quay.io/biocontainers/igv-reports:1.16.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(sites), path(tracks), path(tracks_indices)
    tuple val(meta2), path(fasta), path(fai)
    path  ideogram, stageAs: 'ideogram_in/*'  // optional locally-staged cytoband file; [] to fetch from UCSC at runtime
    path  gnomad_sv_track, stageAs: 'gnomad_tracks/*'  // optional population-SV annotation track (bgzipped VCF or BED); [] to skip
    path  gnomad_sv_tbi,   stageAs: 'gnomad_tracks/*'  // tabix index for gnomad_sv_track; [] when not provided

    output:
    tuple val(meta), path("*.html") , emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Reference mode selector (local copy of the nf-core IGVREPORTS module).
    //
    //   - When `task.ext.genome` is set (e.g. 'hg38'), `create_report` is
    //     invoked with `--genome <id>` which pulls the IGV-hosted genome
    //     refs: RefSeq/GENCODE gene track, chromosome ideogram, etc
    //   - When `task.ext.genome` is unset, use the `--fasta <file>`
    //     report has no gene track and no ideogram unless they are also
    //     staged via `--tracks` / `--ideogram` in `ext.args`.
    //
    // Track-config selector:
    //   When `task.ext.color_bam_by_hp` is true, the module builds
    //   a `tracks_config.json` and passes it via `--track-config`. This
    //   gives per-track properties, HP-tag colorBy/groupBy on alignment, etc
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def genome_id = task.ext.genome ?: ''
    def color_bam_by_hp = task.ext.color_bam_by_hp != null ? task.ext.color_bam_by_hp : false
    def ref_opt   = genome_id ? "--genome ${genome_id}" : (fasta ? "--fasta ${fasta}" : '')
    // When `ideogram` is staged, pass it via --ideogram to bypass the URL fetch
    def ideogram_opt = ideogram ? "--ideogram ${ideogram}" : ''

    // Partition tracks by file extension so we can give each type the right IGV track-config properties. 
    def bams = tracks ? tracks.findAll {
        it.toString().toLowerCase().endsWith('.bam')
    } : []
    def beds = tracks ? tracks.findAll {
        def n = it.toString().toLowerCase()
        n.endsWith('.bed') || n.endsWith('.bed.gz')
    } : []
    def others = tracks ? tracks.findAll {
        def n = it.toString().toLowerCase()
        !(n.endsWith('.bam') || n.endsWith('.bed') || n.endsWith('.bed.gz'))
    } : []

    // BAM coloring: tag-based colorBy + groupBy on the HP haplotag.
    def bam_color_kv = color_bam_by_hp
        ? ',"colorBy":"tag","colorByTag":"HP","groupBy":"tag","groupByTag":"HP"'
        : ''

    // Long-read SV rendering options for alignment tracks (see conf/modules.config
    // IGVREPORTS block for rationale). Each is opt-in via task.ext; defaults are
    // sensible for short-read / SNV reports (no change from prior behaviour).
    def bam_show_soft_clips = task.ext.bam_show_soft_clips != null ? task.ext.bam_show_soft_clips : false
    def bam_display_mode    = task.ext.bam_display_mode    ?: ''
    def bam_sampling_depth  = task.ext.bam_sampling_depth  ?: 0
    def bam_extra_kv = ''
    if (bam_show_soft_clips)          bam_extra_kv += ',"showSoftClips":true'
    if (bam_display_mode)             bam_extra_kv += ',"displayMode":"' + bam_display_mode + '"'
    if ((bam_sampling_depth as int) > 0) bam_extra_kv += ',"samplingDepth":' + bam_sampling_depth

    def entries = []
    bams.each {
        entries << '{"name":"' + it.name + '","url":"' + it.name + '","type":"alignment","format":"bam"' + bam_color_kv + bam_extra_kv + '}'
    }
    // Explicit RefSeq gene track between BAM and BED layers. The genome
    // bundle's default gene track typically renders at the bottom of user tracks
    if (genome_id == "hg38") {
        entries << '{"name":"RefSeq Genes (hg38)","url":"https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.txt.gz","indexURL":"https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.txt.gz.tbi","type":"annotation","format":"refgene","displayMode":"EXPANDED"}'
    } else if (genome_id == "hg19") {
        entries << '{"name":"RefSeq Genes (hg19)","url":"https://s3.amazonaws.com/igv.org.genomes/hg19/refGene.sorted.txt.gz","indexURL":"https://s3.amazonaws.com/igv.org.genomes/hg19/refGene.sorted.txt.gz.tbi","type":"annotation","format":"refgene","displayMode":"EXPANDED"}'
    }
    beds.each {
        entries << '{"name":"' + it.name + '","url":"' + it.name + '","type":"annotation","format":"bed"}'
    }
    others.each {
        entries << '{"name":"' + it.name + '","url":"' + it.name + '"}'
    }

    // Optional gnomAD SV population-frequency track. Rendered as an annotation ribbon (VCF: per-record markers with hover metadata; BED: interval bars).
    if ( gnomad_sv_track ) {
        def gnomad_name = gnomad_sv_track.name
        def gnomad_format = gnomad_name.toLowerCase().endsWith('.bed.gz')
            ? 'bed'
            : (gnomad_name.toLowerCase().endsWith('.vcf.gz') ? 'vcf' : 'annotation')
        def gnomad_index_url = gnomad_sv_tbi ? gnomad_sv_tbi.name : ''
        def index_kv = gnomad_index_url ? ',"indexURL":"' + gnomad_index_url + '"' : ''
        entries << '{"name":"gnomAD SV v4.1","url":"' + gnomad_name + '"' + index_kv +
                   ',"type":"annotation","format":"' + gnomad_format + '","displayMode":"COLLAPSED"}'
    }

    def use_track_config = entries.size() > 0
    def track_config_arg = use_track_config ? "--track-config tracks_config.json" : ""
    def write_track_config = use_track_config
        ? "echo '[" + entries.join(',') + "]' > tracks_config.json"
        : ":"

    // --tracks fallback when no track-config is generated
    def track_arg = (!use_track_config && tracks)
        ? "--tracks " + tracks.collect { track -> track.toString() }.join(' ')
        : ""
    if (args.contains("--tracks") && track_arg) {
        args = args.replace("--tracks", track_arg)
        track_arg = ""
    }

    """
    ${write_track_config}

    create_report $sites \\
        $args \\
        $ref_opt \\
        $ideogram_opt \\
        $track_config_arg \\
        $track_arg \\
        --output ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """
}
