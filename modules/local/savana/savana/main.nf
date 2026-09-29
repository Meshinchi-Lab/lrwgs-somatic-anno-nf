process SAVANA_SAVANA {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/savana:1.3.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pon_1kg), val(min_support), val(a_min_reads), val(cn_step_ch)
    tuple val(meta2), path(fasta), path(fai)
    path(contigs)

    output:
    tuple val(meta), path("**/*.classified.somatic.vcf"), emit: vcf
    tuple val(meta), path("**/*.classified.somatic.bedpe"), emit: bed
    tuple val(meta), path("**/*.{tsv,vcf,bed}"), emit: results
    tuple val("${task.process}"), val('savana'), eval("savana --version 2>&1 | head -1 || echo 'unknown'"), topic: versions, emit: versions_savana

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MPLCONFIGDIR=\$PWD

    # SAVANA's tumour-only `to` subcommand includes copy-number fitting and with no flag to disable it 
    # CN fitting runs AFTER SV calling, so a fitting rejection such as
    #     "Fit of purity=0.32 and ploidy=3.16 is NOT an acceptable solution,
    #      as main CN change step = 5"
    #     "No fits found. See No_fit_found_PARAMS_out.tsv in output"
    # exits 1 with the SV VCF already written. That failure is spurious for this pipeline and the copy-number products are never read.
    #

    savana_rc=0
    savana to \\
        $args \\
        --threads ${task.cpus} \\
        --cna_threads ${task.cpus} \\
        --length 50 \\
        --sample $prefix \\
        --tumour $bam \\
        --outdir ${prefix}_out \\
        --ref $fasta \\
        --g1000_vcf $pon_1kg \\
        --min_support $min_support \\
        --contigs $contigs \\
        --allele_min_reads $a_min_reads \\
        --main_cn_step_change $cn_step_ch \\
        --overwrite || savana_rc=\$?

    # check the SAVANA outputs and convert the error code from 1 to 0 if CN detection is the source of the failure. 
    # if SV results are empty or don't exist, then use an exit 1 for a real pipeline failure, since SAVANA is a required dependency for the SV annotations
    sv_vcf=\$(find ${prefix}_out -name '*.classified.somatic.vcf' -size +0c 2>/dev/null | head -1)
    if [ -z "\$sv_vcf" ]; then
        echo "ERROR: savana exited \$savana_rc and produced no non-empty *.classified.somatic.vcf" >&2
        if [ "\$savana_rc" -ne 0 ]; then exit "\$savana_rc"; fi
        exit 1
    fi
    if [ "\$savana_rc" -ne 0 ]; then
        echo "WARN: savana exited \$savana_rc — copy-number fitting failed (see No_fit_found_PARAMS*.tsv)." >&2
        echo "WARN: SV calls in \$sv_vcf were written before the CNA step and are retained; CNA output is not used by this pipeline." >&2
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.classified.somatic.vcf
    touch ${prefix}/${prefix}.classified.somatic.bedpe
    """
}
