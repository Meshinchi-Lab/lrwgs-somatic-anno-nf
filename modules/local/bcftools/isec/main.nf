// modules/local/bcftools/isec/main.nf
//
// Local bcftools/isec wrapper for 2-input strict intersection. Emits the
// first-input subset (0000.vcf.gz) and second-input subset (0001.vcf.gz) as
// distinct named output channels so the downstream subworkflow doesn't have
// to crack open the isec output directory and pick files by name.
//
// Inputs:
//   tuple val(meta), path(vcfs), path(tbis)
//     vcfs = [first_input.vcf.gz, second_input.vcf.gz]
//     tbis = [first_input.vcf.gz.tbi, second_input.vcf.gz.tbi]
//
// Outputs:
//   tuple val(meta), path("isec/0000.vcf.gz"), path("isec/0000.vcf.gz.tbi"), emit: a
//     — records from FIRST input shared by both, exact REF+ALT match
//   tuple val(meta), path("isec/0001.vcf.gz"), path("isec/0001.vcf.gz.tbi"), emit: b
//     — records from SECOND input shared by both, exact REF+ALT match
//
// Command flags (set via task.ext.args in modules.config):
//   -n+2          — output positions present in ≥ 2 input files
//   -c none       — exact REF+ALT match (not just position)
//   -Oz           — bgzipped VCF output
//   --write-index=tbi — auto-create .tbi alongside each output
//
process BCFTOOLS_ISEC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("isec/0000.vcf.gz"), path("isec/0000.vcf.gz.tbi"), emit: a
    tuple val(meta), path("isec/0001.vcf.gz"), path("isec/0001.vcf.gz.tbi"), emit: b
    tuple val(meta), path("isec/sites.txt"),  emit: sites,   optional: true
    tuple val(meta), path("isec/README.txt"), emit: readme,  optional: true
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-n+2 -c none -Oz --write-index=tbi'
    """
    bcftools isec ${args} -p isec ${vcfs.join(' ')}
    """

    stub:
    """
    mkdir -p isec
    touch isec/0000.vcf.gz isec/0000.vcf.gz.tbi
    touch isec/0001.vcf.gz isec/0001.vcf.gz.tbi
    touch isec/sites.txt isec/README.txt
    """
}
