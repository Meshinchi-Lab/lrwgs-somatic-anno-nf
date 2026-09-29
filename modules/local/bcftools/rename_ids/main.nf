// modules/local/bcftools/rename_ids/main.nf
//
// Rewrite VCF ID column to simple sequential `${prefix}_N` strings.
//
// Why: downstream merging in ADD_CALLER joins the AnnotSV TSV with a
// caller-lookup TSV on the `ID` column. AnnotSV's `variantconvert` step
// shifts the VCF POS by +1 (1-based normalization) and, for callers that
// embed coordinates inside the ID (e.g. Wakhan `wakhan:GAIN:chr21:13000001-31877979`),
// annotsv ends up rewriting the the embedded coordinates in the ID. 
//
// Match strategy: Exact CHROM+POS is unambiguous because each input record has a unique CHROM+POS key.
//
process BCFTOOLS_RENAME_IDS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)
    val   id_prefix

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | sed '1!d; s/^.*bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: "${task.ext.args}"
    args2 = task.ext.args2 ?: "${task.ext.args2}"
    prefix = task.ext.prefix ?: "${meta.id}.renamed"
    """
    ### bcftools does not annotate 0-based (telomeres) positions, though 0-based positions do meet valid VCF specs. 
    #Option 1: bcftools annotate --set-id '%SVTYPE\\_%CHROM\\_%POS\\_%INFO/END' -Oz -o ${prefix}.vcf.gz ${vcf}
    
    cat > ID_hdr.txt <<EOF
##INFO=<ID=WAKHAN_ID,Number=1,Type=String,Description="Original ID assigned by WAKHAN">
EOF

    sed -E "s/\\t0{1}\\t${id_prefix}/\\t1\\t${id_prefix}/" ${vcf} > ${prefix}.1pos.vcf

    bcftools query -f '%CHROM\\t%POS\\t%ID\\n' ${prefix}.1pos.vcf \\
        | awk -v p='${id_prefix}' 'BEGIN{OFS="\\t"} { print \$1, \$2, p"_"NR, \$3 }' \\
        | bgzip -c > rename.tab.gz

    # add the -0 flag later, when future compatibility for the position does not need to be temporarily converted
    tabix -s1 -b2 -e2 rename.tab.gz
    
    bcftools annotate \\
        -h ID_hdr.txt \\
        -a rename.tab.gz \\
        -c CHROM,POS,ID,INFO/WAKHAN_ID \\
        ${prefix}.1pos.vcf \\
        | sed -E "s/\\t1{1}\\t${id_prefix}/\\t0\\t${id_prefix}/" \\
        | bcftools view -Oz -o ${prefix}.vcf.gz
    
    tabix -p vcf ${prefix}.vcf.gz
    rm ${prefix}.1pos.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.renamed"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
