process MINDA_EXTRACT_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::htslib=1.21 conda-forge::python=3.11"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_bedtools_htslib_python:4936763617c5639d'}"

    input:
    tuple val(meta), path(minda_vcf),
                     path(severus_vcf), path(severus_tbi),
                     path(savana_vcf),  path(savana_tbi)

    output:
    tuple val(meta), path("*.annotation.tab.gz"), path("*.annotation.tab.gz.tbi"), emit: annotation_vcf
    tuple val(meta), path("*.header.txt"),                                          emit: header
    tuple val(meta), path("*.savana_supported_ids.txt"),                            emit: supported_ids
    tuple val("${task.process}"), val('minda_extract_annotations'), val('1.0.0'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minda_extract_annotations.py \\
        --minda_vcf   ${minda_vcf} \\
        --severus_vcf ${severus_vcf} \\
        --savana_vcf  ${savana_vcf} \\
        --sample_id   ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.annotation.tab.gz
    touch ${prefix}.annotation.tab.gz.tbi
    touch ${prefix}.header.txt
    touch ${prefix}.savana_supported_ids.txt
    """
}
