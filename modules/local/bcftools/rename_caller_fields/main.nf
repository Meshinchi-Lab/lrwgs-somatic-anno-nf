// modules/local/bcftools/rename_caller_fields/main.nf
//
// Renames FORMAT/AF and FORMAT/DP to caller-tagged names so the cross-caller
// VAF and depth can be carried in INFO without name collisions across callers.
//
// Inputs:
//   tuple val(meta), path(vcf), path(tbi)
//   val   caller_prefix    — e.g. "CLAIRSTO" or "DEEPSOMATIC" — produces
//                            FORMAT/{prefix}_VAF and FORMAT/{prefix}_DP
//
// Outputs:
//   tuple val(meta), path("*.renamed.vcf.gz"), path("*.renamed.vcf.gz.tbi")
//
// Behaviour: writes a rename map at runtime; runs `bcftools annotate
// --rename-annots`. The map is FORMAT-level not INFO-level.
//
process BCFTOOLS_RENAME_CALLER_FIELDS {
    tag "${meta.id}/${caller_prefix}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   caller_prefix
    val   source_vaf_field   // "AF" for ClairS-TO, "VAF" for DeepSomatic

    output:
    tuple val(meta), path("*.renamed.vcf.gz"), path("*.renamed.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${caller_prefix.toLowerCase()}"
    """
    # Build rename map: FORMAT/<source_vaf_field> -> FORMAT/<CALLER>_VAF, FORMAT/DP -> FORMAT/<CALLER>_DP
    cat > rename.txt <<EOF
FORMAT/${source_vaf_field} ${caller_prefix}_VAF
FORMAT/DP ${caller_prefix}_DP
EOF

    bcftools annotate \\
        --rename-annots rename.txt \\
        ${vcf} \\
        -Oz -o ${prefix}.renamed.vcf.gz
    bcftools index -t ${prefix}.renamed.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${caller_prefix.toLowerCase()}"
    """
    touch ${prefix}.renamed.vcf.gz
    touch ${prefix}.renamed.vcf.gz.tbi
    """
}
