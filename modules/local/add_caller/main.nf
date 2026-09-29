process ADD_CALLER {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bioframe=0.8.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.8.0--pyhdfd78af_0' :
        'quay.io/biocontainers/bioframe:0.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(annotsv_tsv), path(caller_tsv)
    path candidate_genes   // optional one-gene-per-line file (e.g. OncoKB cancer-gene list); pass [] to omit

    output:
    // Full merged TSV — every split-mode gene row retained. Feeds KNOTANNOTSV_XL.
    tuple val(meta), path("${prefix}.tsv"),       emit: tsv
    // HTML feed — split rows per variant capped at `task.ext.max_split_rows`
    tuple val(meta), path("${prefix}.html.tsv"),  emit: tsv_html
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.slim.annotated"
    def cand_arg = candidate_genes ? "--candidate-genes ${candidate_genes}" : ""
    def max_split_rows = task.ext.max_split_rows ?: 10
    """
    add_caller.py \\
        --annotsv-tsv ${annotsv_tsv} \\
        --caller-tsv  ${caller_tsv} \\
        ${cand_arg} \\
        --output         ${prefix}.tsv \\
        --output-html    ${prefix}.html.tsv \\
        --max-split-rows ${max_split_rows}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.slim.annotated"
    """
    touch ${prefix}.tsv
    touch ${prefix}.html.tsv
    """
}
