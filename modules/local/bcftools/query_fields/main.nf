// modules/local/bcftools/query_fields/main.nf
//
// Generic `bcftools query -f <FORMAT>` wrapper that emits a TSV of extracted
// per-record fields. The FORMAT string is supplied as a val input so aliased
// invocations (BCFTOOLS_QUERY_SV_FIELDS / BCFTOOLS_QUERY_SNV_FIELDS in
// IGV_REPORTS_SV) can pass different projections without touching modules.config.
//
// Used by IGV_REPORTS_SV to feed the Python `bin/vcf_to_bedpe.py` converter
// (which runs in a Python-only container that lacks bcftools).
//
// Inputs:
//   tuple val(meta), path(vcf), path(tbi)
//   val   format_string   — bcftools -f expression, e.g. '%CHROM\\t%POS\\n'
//
// Outputs:
//   tuple val(meta), path("*.tsv"), emit: tsv
//
process BCFTOOLS_QUERY_FIELDS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   format_string

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools query -f '${format_string}' ${vcf} > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
