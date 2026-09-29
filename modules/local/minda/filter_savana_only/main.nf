process MINDA_FILTER_SAVANA_ONLY {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::htslib=1.21 conda-forge::python=3.11"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_bedtools_htslib_python:4936763617c5639d'}"

    input:
    tuple val(meta), path(savana_vcf), path(savana_tbi), path(supported_ids), path(header_txt)

    output:
    tuple val(meta), path("*.savana_only.vcf"), emit: vcf
    tuple val("${task.process}"), val('minda_filter_savana_only'), val('1.0.0'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minda_filter_savana_only.py \\
        --savana_vcf    ${savana_vcf} \\
        --supported_ids ${supported_ids} \\
        --header_txt    ${header_txt} \\
        --sample_id     ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.savana_only.vcf
    """
}
