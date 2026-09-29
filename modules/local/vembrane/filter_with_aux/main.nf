// modules/local/vembrane/filter_with_aux/main.nf
//
// Thin wrapper over vembrane filter that adds an explicit `path(aux_file)`
// input. The nf-core vembrane/filter module doesn't have an aux file input,
// so we duplicate its structure here and pass --aux NAME=PATH alongside the
// expression. Used by FILTER_SNV to make the candidate gene list available
// as AUX["candidate_genes"] in the filter expression.
//
// Inputs:
//   tuple val(meta), path(vcf), path(tbi)
//   val   expression       — the vembrane filter expression
//   path  aux_file         — file with one entry per line (e.g. gene symbols).
//                            Pass [] to skip (and the --aux flag is omitted).
//   val   aux_name         — the AUX key name, e.g. 'candidate_genes'
//
// Outputs:
//   tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi")
//
process VEMBRANE_FILTER_WITH_AUX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vembrane:2.4.0--pyhdfd78af_0' :
        'quay.io/biocontainers/vembrane:2.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   expression
    path  aux_file
    val   aux_name

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('vembrane'),
          eval("vembrane --version | sed '1!d;s/.* //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}.snv"
    def aux_arg = aux_file ? "--aux ${aux_name}=${aux_file}" : ''
    """
    # vembrane emits an uncompressed VCF; the downstream BCFTOOLS_SORT step in FILTER_SNV handles bgzip + tabix indexing via `-Oz --write-index=tbi`.
    vembrane filter \\
        ${aux_arg} \\
        ${args} \\
        ${expression} \\
        ${vcf} \\
        -o ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.snv"
    """
    touch ${prefix}.vcf
    """
}
