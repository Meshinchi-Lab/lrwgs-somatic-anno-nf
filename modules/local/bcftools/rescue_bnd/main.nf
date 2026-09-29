process BCFTOOLS_RESCUE_BND {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    // filtered_vcf: uncompressed VCF output from vembrane (may have orphaned BND partners)
    // original_vcf: the VCF that was passed to vembrane (VEP-annotated or AnnotSV-annotated)
    tuple val(meta), path(filtered_vcf), path(original_vcf)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -1 | sed 's/bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.merged.annotated.rescued"

    """
    # Collect IDs of every record that passed vembrane
    # bcftools view ${filtered_vcf}
    bcftools query -f '%ID\\n' ${filtered_vcf} | sort -u > passing_ids.txt

    # From the original VCF, find BND records whose MATEID partner passed the
    # filter but whose own ID did not — these are orphaned breakend partners.
    # awk logic:
    #   first file  (passing_ids.txt): build pass[] lookup
    #   second input (bcftools query): col1=ID, col2=MATE_ID
    #     keep rows where MATE_ID is in pass[] but the record ID itself is not
    bcftools view -i 'SVTYPE="BND"' ${original_vcf} \\
        | bcftools query -f '%ID\\t%INFO/MATE_ID\\n' \\
        | awk 'NR==FNR{pass[\$1]=1; next} \$2!="." && (\$2 in pass) && !(\$1 in pass)' \\
              passing_ids.txt - \\
        | cut -f1 | sort -u > rescue_ids.txt

    if [ -s rescue_ids.txt ]; then
        bcftools view -i 'ID=@rescue_ids.txt' ${original_vcf} -o rescued_partners.vcf
        bgzip ${filtered_vcf} && tabix ${filtered_vcf}.gz
        bgzip rescued_partners.vcf && tabix rescued_partners.vcf.gz
        bcftools concat -a ${filtered_vcf}.gz rescued_partners.vcf.gz | bcftools sort -o ${prefix}.vcf
    else
        cp ${filtered_vcf} ${prefix}.vcf
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.merged.annotated.rescued"
    """
    touch ${prefix}.vcf
    """
}
