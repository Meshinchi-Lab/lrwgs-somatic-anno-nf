// modules/local/sv_to_bedpe/main.nf
//
// consume the pre-extracted TSV fields from
// BCFTOOLS_QUERY_SV_FIELDS (and optionally BCFTOOLS_QUERY_SNV_FIELDS) and emit
// three bucket-split BEDPE files. 
//
// The Python script (bin/vcf_to_bedpe.py) does all string manipulation —
// BND ALT bracket parsing, VEP CSQ SYMBOL extraction, INS SVLEN +Nbp name substitution, inter-anchor-distance bucket routing.
//
// Inputs:
//   tuple val(meta), path(sv_tsv), path(snv_tsv)
//     - snv_tsv is `[]` sentinel when SNVs are disabled; script skips --snv-in.
//
// Outputs:
//   tuple val(meta), path("*.pass2.sites.small.bedpe"),  emit: small
//   tuple val(meta), path("*.pass2.sites.medium.bedpe"), emit: medium
//   tuple val(meta), path("*.pass2.sites.large.bedpe"),  emit: large
//
process SV_TO_BEDPE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    // Full python:3.11 image (not -slim) — includes `procps` which Nextflow requires for task metrics collection (`ps` command). 
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.11' :
        'docker.io/python:3.11' }"

    input:
    tuple val(meta), path(sv_tsv), path(snv_tsv, stageAs: 'snv_input/*')
    // Optional one-gene-per-line candidate-gene panel. `[]` sentinel from the
    // subworkflow when params.candidateGenesFile is not set — script skips
    // --candidate-genes cleanly and falls back to first-10-of-Gene_name.
    path candidate_genes, stageAs: 'candidate_genes/*'

    output:
    tuple val(meta), path("*.pass2.sites.small.bedpe"),  emit: small
    tuple val(meta), path("*.pass2.sites.medium.bedpe"), emit: medium
    tuple val(meta), path("*.pass2.sites.large.bedpe"),  emit: large
    tuple val("${task.process}"), val('python'),
          eval("python3 --version | sed 's/Python //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snv_flag = snv_tsv ? "--snv-in ${snv_tsv}" : ""
    def cand_flag = candidate_genes ? "--candidate-genes ${candidate_genes}" : ""
    """
    python3 ${projectDir}/bin/vcf_to_bedpe.py \\
        --sv-in ${sv_tsv} ${snv_flag} ${cand_flag} \\
        --out-small  ${prefix}.pass2.sites.small.bedpe \\
        --out-medium ${prefix}.pass2.sites.medium.bedpe \\
        --out-large  ${prefix}.pass2.sites.large.bedpe
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pass2.sites.small.bedpe
    touch ${prefix}.pass2.sites.medium.bedpe
    touch ${prefix}.pass2.sites.large.bedpe
    """
}
