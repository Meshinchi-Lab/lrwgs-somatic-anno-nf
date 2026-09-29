// modules/local/bcftools/add_caller_tag/main.nf
//
// Injects an INFO/CALLER tag with a static value across every record in a VCF.
// Used by MERGE_SNV to tag the consensus VCF with "CLAIRSTO,DEEPSOMATIC" so
// downstream filters / reports can trace provenance.
//
// Inputs:
//   tuple val(meta), path(vcf), path(tbi)
//   val   caller_value   — e.g. "CLAIRSTO,DEEPSOMATIC"
//
// Outputs:
//   tuple val(meta), path("*.caller.vcf.gz"), path("*.caller.vcf.gz.tbi")
//
process BCFTOOLS_ADD_CALLER_TAG {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   caller_value
    val   id_prefix

    output:
    tuple val(meta), path("*.caller.vcf.gz"), path("*.caller.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Inject the INFO/CALLER header line if not already present, then set the
    # value on every record using a synthetic per-position annotation source.
    cat > caller_hdr.txt <<EOF
##INFO=<ID=CALLER,Number=1,Type=String,Description="Caller(s) that emitted this record">
EOF
    if [[ "${meta.variant_type}" == "cnv" ]] 
    then
        ### bcftools does not annotate 0-based (telomeres) positions, though 0-based positions do meet valid VCF specs. 
       bcftools view ${vcf} | sed -E "s/\\t0{1}\\t${id_prefix}/\\t1\\t${id_prefix}/"  > ${prefix}.1pos.vcf && mv ${prefix}.1pos.vcf ${vcf}
    fi

    # Build an annotation source: a TAB-separated file with CHROM POS REF ALT VALUE.
    bcftools query \\
        -f '%CHROM\\t%POS\\t${caller_value}\\n' \\
        ${vcf} \\
        | bgzip -c > annot.tsv.gz
    tabix -s1 -b2 -e2 annot.tsv.gz

    bcftools annotate \\
        --header-lines caller_hdr.txt \\
        -a annot.tsv.gz \\
        -c CHROM,POS,INFO/CALLER \\
        ${vcf} > file.vcf
    
    if [[ "${meta.variant_type}" == "cnv" ]] 
    then
        ### bcftools does not annotate 0-based (telomeres) positions, though 0-based positions do meet valid VCF specs. 
       sed -E "s/\\t0{1}\\t${id_prefix}/\\t1\\t${id_prefix}/" file.vcf > file.0pos.vcf && mv file.0pos.vcf file.vcf
    fi
    
    # convert to gzipped format and index the file
    bcftools view -Oz -o ${prefix}.caller.vcf.gz file.vcf
    bcftools index -t ${prefix}.caller.vcf.gz

    rm annot.tsv.gz annot.tsv.gz.tbi caller_hdr.txt file*vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.caller.vcf.gz
    touch ${prefix}.caller.vcf.gz.tbi
    """
}

/*
// modules/local/bcftools/add_caller_tag/main.nf
//
// Inject a constant INFO/CALLER tag into every record of a VCF.
// The per-record annotation table is self-derived from the input VCF via
// `bcftools query`, then `bcftools annotate` applies it. Pre-verified at the
// CLI (see plan task 1) before being wrapped here.
//
process BCFTOOLS_ADD_CALLER_TAG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)
    val   caller_name

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | sed '1!d; s/^.*bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.caller_tag"
    """
    echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller (constant)">' > caller.hdr

    bcftools query -f '%CHROM\\t%POS\\t${caller_name}\\n' ${vcf} \\
        | bgzip -c > caller.tab.gz
    tabix -s1 -b2 -e2 caller.tab.gz

    bcftools annotate \\
        -a caller.tab.gz \\
        -h caller.hdr \\
        -c CHROM,POS,CALLER \\
        -Oz \\
        -o ${prefix}.vcf.gz \\
        ${vcf}
    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.caller_tag"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
*/
