// Remove from a VCF every record whose position is already present in another VCF.
//
// Used by the candidate-gene rescue in MERGE_SNV: a rescued single-caller record
// must NOT be concatenated alongside the consensus record for the same position.
//
// Why this exists rather than relying on `bcftools concat -D`: with `-a -D`,
// concat merges records in position order across all inputs and drops duplicates
// by whichever copy it encounters first. That is NOT guaranteed to be the
// consensus copy. In practice it kept the rescue copy for 2,486 positions, which
// silently downgraded genuine two-caller consensus calls to single-caller
// (CALLER=CLAIRSTO,DEEPSOMATIC -> CALLER=CLAIRSTO) and discarded one caller's
// VAF/DP. Subtracting first makes `-D` belt-and-braces instead of load-bearing.
//
// `-T ^FILE` is the negated targets form: keep records NOT in FILE.
process BCFTOOLS_SUBTRACT_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(exclude), path(exclude_tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}.subtracted"
    def args = task.ext.args ?: ''
    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    bcftools view ${args} -T ^${exclude} ${vcf} -Oz -o ${prefix}.vcf.gz
    bcftools index -t ${prefix}.vcf.gz

    echo "INFO: \$(bcftools index -n ${vcf}) records in, \$(bcftools index -n ${prefix}.vcf.gz) kept after removing positions present in ${exclude}" >&2
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.subtracted"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
