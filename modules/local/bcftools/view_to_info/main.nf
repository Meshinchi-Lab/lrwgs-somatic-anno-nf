// modules/local/bcftools/view_to_info/main.nf
//
// Promotes a list of FORMAT fields to INFO using `bcftools annotate
// -c INFO/X:=FORMAT/X`. The colon-equals operator copies the per-sample
// FORMAT value to the top-level INFO column. Since our input is single-sample,
// this is a clean lift; multi-sample VCFs would need a per-sample annotation
// strategy instead (out of scope for this pipeline).
//
// Inputs:
//   tuple val(meta), path(vcf), path(tbi)
//   val   field_names   — list of FORMAT field names to promote (e.g.
//                         ['CLAIRSTO_VAF', 'CLAIRSTO_DP']). They are joined
//                         into the -c argument as INFO/X:=FORMAT/X pairs.
//
// Outputs:
//   tuple val(meta), path("*.info.vcf.gz"), path("*.info.vcf.gz.tbi")
//
process BCFTOOLS_VIEW_TO_INFO {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   field_names

    output:
    tuple val(meta), path("*.info.vcf.gz"), path("*.info.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Build the bcftools query -f format string: %CHROM\t%POS\t%REF\t%ALT[\t%FIELD1\t%FIELD2...]
    def query_fields = field_names.collect { "%${it}" }.join('\\t')
    // Build the annotate -c column spec
    def info_columns = field_names.collect { "INFO/${it}" }.join(',')
    """
    # Step 1: Extract per-sample FORMAT values to a tab-separated annotation source.
    # The bracket [] notation extracts per-sample FORMAT fields; for single-sample
    # VCFs (our case) this is one value per field per record.
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t${query_fields}]\\n' \\
        ${vcf} \\
        | bgzip -c > annot.tsv.gz
    tabix -s1 -b2 -e2 annot.tsv.gz

    # Step 2: Synthesize INFO header lines from the source FORMAT header lines.
    bcftools view -h ${vcf} | grep -E "^##FORMAT=<ID=(${field_names.join('|')})," \\
        | sed 's/^##FORMAT=/##INFO=/' \\
        > new_info_headers.txt

    # Step 3: Annotate the original VCF with the extracted values, now at INFO level.
    bcftools annotate \\
        --header-lines new_info_headers.txt \\
        -a annot.tsv.gz \\
        -c CHROM,POS,REF,ALT,${info_columns} \\
        ${vcf} \\
        -Oz -o ${prefix}.info.vcf.gz
    bcftools index -t ${prefix}.info.vcf.gz

    rm annot.tsv.gz annot.tsv.gz.tbi new_info_headers.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.info.vcf.gz
    touch ${prefix}.info.vcf.gz.tbi
    """
}
