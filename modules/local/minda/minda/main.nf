process MINDA_MINDA {
    tag "$meta.id"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    container 'quay.io/jennylsmith/minda:1.0'

    input:
    tuple val(meta), path(vcf1), path(vcf2)
    val(tolerance)

    output:
    //example tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("**/*ensemble.vcf"), emit: ensemble_vcf
    tuple val(meta), path("**/*.tsv"), emit: metrics
    tuple val(meta), path("**/*.log"), emit: log
    tuple val(meta), path("**/results/*.tsv"), emit: summaries
    
    // tuple val("${task.process}"), val('minda'), eval("minda --version"), topic: versions, emit: versions_minda

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   """
    minda ensemble \\
        --sample_name ${prefix} \\
        --vcfs $vcf1 $vcf2 \\
        --out_dir ${prefix} \\
        --min_support 2 \\
        --tolerance ${tolerance}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/results
    touch ${prefix}/${prefix}_minda_ensemble.vcf
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.log
    touch ${prefix}/results/${prefix}.tsv
    """
}
