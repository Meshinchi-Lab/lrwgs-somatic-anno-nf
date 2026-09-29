// Make a VCF acceptable to strict htsjdk-based tools (picard, GATK).
//
// htslib/bcftools tolerate two things the VCF spec forbids and htsjdk rejects
// outright. Both occur in the CIViC reference VCF:
//
//   1. Whitespace inside INFO values. VCF 4.2 §1.4.1 forbids whitespace,
//      semicolons and equals-signs in INFO. CIViC embeds free text such as
//      "Missense Variant", "Skin Melanoma" and "Vemurafenib (NCIt ID C64768)".
//      picard: "The VCF specification does not allow for whitespace in the
//      INFO field".
//   2. Records where REF == ALT (CIViC has one: chr7:148524752 C>C, the EZH2
//      "Intron 6 Mutation" region entry). picard: "Duplicate allele added to
//      VariantContext". Such a record can never match a real variant call, so
//      dropping it loses nothing.
//
// Only (2) has a native bcftools expression (`-e 'REF=ALT'`, supplied via
// ext.args). bcftools has no value-level string substitution — annotate can
// transfer, remove and rename tags, but cannot rewrite characters inside a
// value — so (1) needs one `sed` restricted to non-header lines. That is safe
// and total: VCF columns are TAB-separated and every data column forbids
// spaces, so a space on a data line is always illegal and always substitutable.
// Header lines are untouched, preserving the human-readable INFO Descriptions.
process BCFTOOLS_SANITIZE_VCF {
    tag "${meta.id ?: 'vcf'}"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id ?: 'input'}.sanitized"
    def args = task.ext.args ?: ''
    if ("${vcf}" == "${prefix}.vcf.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # `-e 'REF=ALT'` (ext.args) drops spec-invalid duplicate-allele records natively.
    # sed then replaces every space on a data line with "_"; header lines keep theirs.
    bcftools view ${args} ${vcf} \\
        | sed '/^#/!s/ /_/g' \\
        | bcftools view -Oz -o ${prefix}.vcf.gz

    bcftools index -t ${prefix}.vcf.gz

    # Fail loudly rather than hand a still-invalid file to picard downstream.
    if bcftools view -H ${prefix}.vcf.gz | grep -q ' '; then
        echo "ERROR: whitespace remains on a data line after sanitization." >&2
        exit 1
    fi
    echo "INFO: sanitized \$(bcftools index -n ${prefix}.vcf.gz) records into ${prefix}.vcf.gz"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id ?: 'input'}.sanitized"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
