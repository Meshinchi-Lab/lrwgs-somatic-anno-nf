process STRIP_FLAG_CSQ {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::htslib=1.21"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.21--h566b1c6_1' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.csq_stripped"
    """
    strip_flag_csq.sh ${vcf} ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.csq_stripped"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    """
}
