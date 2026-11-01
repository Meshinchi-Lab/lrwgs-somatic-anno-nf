// modules/local/circos/main.nf
// Render ONE circos figure (PNG + SVG) from a fixed-name track set with the
// repo's assets/genomic_viz/circos configuration.
//
// Inputs are staged under FIXED names because tracks.conf references
// links.txt / cna.txt / labels.txt literally. The configuration is passed as the
// five .conf FILES (staged under conf/), not as the directory: Nextflow
// fingerprints a directory input by path, size and mtime only, so `-resume`
// would silently reuse a figure after a .conf edit. The hg38 karyotype ships
// inside the circos container and circos.conf resolves it through circos's own
// data search path, so nothing install-specific is needed here.
//
// A zero-byte track file switches that track off through `-param show_<x>=no`
// (tracks.conf reads `show = conf(show_<x>)`). circos 0.69 renders empty data
// files without error anyway, so the switch is defensive and keeps the log
// explicit about which tracks a figure carries.
//
// Contract: meta.layers (positive integer, default 1) is the tile stack depth --
// 1 for a per-sample figure, 6 for the cohort ring. It is validated here and left
// in the emitted meta; consumers collect the outputs by file name, not by meta.
process CIRCOS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::circos=0.69.9"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circos:0.69.9--hdfd78af_0' :
        'quay.io/biocontainers/circos:0.69.9--hdfd78af_0' }"

    input:
    tuple val(meta), path('links.txt'), path('cna.txt'), path('labels.txt')
    path conf, stageAs: 'conf/*'    // the assets/genomic_viz/circos/*.conf FILES (a list)

    output:
    tuple val(meta), path("${prefix}.png"), emit: png
    tuple val(meta), path("${prefix}.svg"), emit: svg
    // bioconda circos 0.69.9 still reports itself as "0.69-8"; this records the
    // binary's own answer rather than the package tag.
    tuple val("${task.process}"), val('circos'),
          eval("circos -v | sed -E 's/^circos \\| v ([^ ]+) .*/\\1/'"),
          topic: versions, emit: versions_circos

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.circos"
    def layers_raw = meta.layers == null ? '1' : meta.layers.toString()
    if ( !layers_raw.isInteger() || layers_raw.toInteger() < 1 ) {
        error("CIRCOS: meta.layers must be a positive integer, got '${meta.layers}' for ${meta.id}")
    }
    def layers = layers_raw.toInteger()
    """
    # Defensive track switches (see header).
    show_sv_links=yes; [ -s links.txt  ] || show_sv_links=no
    show_cna=yes;      [ -s cna.txt    ] || show_cna=no
    show_labels=yes;   [ -s labels.txt ] || show_labels=no
    echo "INFO: tracks sv_links=\${show_sv_links} cna=\${show_cna} labels=\${show_labels} tile_layers=${layers}" >&2

    # Drop only circos's known harmless Perl warning. The subshell's `|| true` stops
    # grep -v from failing the pipeline when it filters every line; pipefail
    # (process.shell in nextflow.config) still propagates circos's own exit status,
    # so an out-of-memory 137 reaches the retry policy in base.config.
    circos \\
        -conf conf/circos.conf \\
        -param show_sv_links=\${show_sv_links} \\
        -param show_cna=\${show_cna} \\
        -param show_labels=\${show_labels} \\
        -param tile_layers=${layers} \\
        -outputdir . \\
        -outputfile ${prefix}.png \\
        -silent 2>&1 | { grep -v 'Configuration.pm line 781' >&2 || true; }

    # circos streams the SVG while drawing and writes the PNG last: a crash leaves
    # a partial SVG and no PNG, so a non-empty PNG marks a completed run.
    [ -s ${prefix}.png ] || { echo "ERROR: circos wrote no ${prefix}.png" >&2; exit 1; }
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.circos"
    """
    printf 'stub' > ${prefix}.png
    printf '<svg xmlns="http://www.w3.org/2000/svg"></svg>' > ${prefix}.svg
    """
}
