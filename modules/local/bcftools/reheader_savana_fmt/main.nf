// Add missing SAVANA FORMAT header definitions to the MINDA merged VCF.
// MINDA does not carry over per-caller FORMAT ##FORMAT lines, so VAF, hVAF, DR, and DV are present in records but undefined in the header, causing
// bcftools to warn and assume Type=String. bcftools annotate --header-lines is additive — it appends only, without replacing the existing header.
process BCFTOOLS_REHEADER_SAVANA_FMT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.ensemble"

    """
    # hVAF is declared Number=3,Type=Float (HP1, HP2, unphased) — SEVERUS emits three comma-separated values per record
    cat > savana_fmt_headers.txt << 'HDRS'
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Tumour variant allele frequency">
##FORMAT=<ID=hVAF,Number=3,Type=Float,Description="Haplotype-aware tumour variant allele frequency (HP1,HP2,unphased ONT haplotag)">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Tumour reference-supporting reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Tumour variant-supporting reads">
HDRS

    bcftools annotate \\
        --header-lines savana_fmt_headers.txt \\
        -Oz \\
        -o ${prefix}.vcf.gz \\
        ${vcf}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.ensemble"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    """
}
