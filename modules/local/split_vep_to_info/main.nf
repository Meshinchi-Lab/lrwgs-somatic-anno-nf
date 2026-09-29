// modules/local/split_vep_to_info/main.nf
//
// Extract selected VEP CSQ subfields into top-level INFO/* tags and drop
// INFO/CSQ entirely — bypasses vembrane's CSQ parser to sidestep the
// documented round-trip corruption on very long/complex CSQ values
// (notes.2026.08.06.md §3.1, ID_10322_1 fixture). See spec 2026-07-07-vep-csq-bypass-design.md.
//
// Mechanism: `bcftools +split-vep -c <fields>` promotes each named CSQ
// subfield to a top-level INFO tag with the declared Type. Piped into
// `bcftools annotate -x INFO/CSQ` which strips both the header line and
// every record's CSQ value, so downstream vembrane sees no CSQ at all.
//
// Fields to extract are configured per-invocation via `task.ext.split_vep_fields`,
// e.g. 'gnomAD_SV:String,gnomAD_SV_AF:Float' (default in modules.config).
//
process SPLIT_VEP_TO_INFO {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.split_vep.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.split_vep.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fields = task.ext.split_vep_fields ?: 'gnomAD_SV:String,gnomAD_SV_AF:Float'
    """
    # Two-stage pipe: +split-vep extracts CSQ subfields into INFO/*,
    # annotate -x removes INFO/CSQ header + values. Both stages are stream-safe.
    bcftools +split-vep -c '${fields}' -O v ${vcf} \\
        | bcftools annotate -x INFO/CSQ -O z -o ${prefix}.split_vep.vcf.gz
    tabix -p vcf ${prefix}.split_vep.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.split_vep.vcf.gz
    touch ${prefix}.split_vep.vcf.gz.tbi
    """
}
