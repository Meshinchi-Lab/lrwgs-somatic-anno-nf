// Build rescue intervals from a list of gene symbols.
//
// Source is the AnnotSV RefSeq gene BED already staged by PREPARE_REFERENCES:
//   col1 chrom (NCBI-style, unprefixed), col2 start, col3 end, col4 strand, col5 symbol
//
// Three transforms, all required:
//   1. filter to the requested symbols
//   2. strand-aware upstream padding -- lower coordinate for '+', HIGHER for '-'.
//      Not optional: TERT is on the minus strand with a gene body ending at
//      5:1295068, while the C228T promoter hotspot is at chr5:1295113, 45 bp
//      past the higher coordinate. A non-strand-aware pad extends the wrong end.
//   3. chr-prefix, sort, merge overlaps -- the BED ships '5', the callset uses 'chr5'
process BUILD_RESCUE_BED {
    tag "${genes.simpleName}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity','apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(annotations)
    path genes
    val  genome_build
    val  pad

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "${prefix}.dropped.txt",           emit: dropped
    tuple val("${task.process}"), val('bedtools'),
          eval("bedtools --version | sed 's/bedtools //'"),
          topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: 'rescue_regions'
    """
    GENE_BED="${annotations}/Annotations_Human/Genes/${genome_build}/genes.RefSeq.sorted.bed"
    if [ ! -f "\$GENE_BED" ]; then
        echo "ERROR: gene BED not found at \$GENE_BED" >&2
        exit 1
    fi

    # Non-gene entries in the ST9 list (e.g. 'CNVgains 8,10,19', 'Xq', 'Xp22.33')
    # never match column 5; record them so the count is auditable.
    #
    # Many dropped entries are historical gene ALIASES (A20->TNFAIP3, ABL->ABL1,
    # ABCB2->TAP1, ...). The RefSeq BED uses current HGNC symbols only, so aliases
    # never match column 5 -- but the canonical symbol is also present in the list
    # and does match, so no gene is actually lost. Verified for all T-ALL drivers
    # (TERT, BRAF, NOTCH1, CDKN2A, PTEN, FBXW7, PHF6, WT1, RUNX1, IL7R, JAK3,
    # STIL, TAL1, LMO2, TLX3). Expect ~3400 of 4803 entries to drop on the ST9 list.
    awk -F'\\t' 'NR==FNR{g[\$1];next} \$5 in g {seen[\$5]=1} END{for(k in g) if(!(k in seen)) print k}' \\
        ${genes} "\$GENE_BED" | sort > ${prefix}.dropped.txt

    awk -F'\\t' -v pad=${pad} 'NR==FNR{g[\$1];next}
        \$5 in g {
            s=\$2; e=\$3;
            if (\$4=="-") e=e+pad; else s=(s>pad ? s-pad : 0);
            print "chr"\$1"\\t"s"\\t"e"\\t"\$5
        }' ${genes} "\$GENE_BED" \\
      | sort -k1,1 -k2,2n \\
      | bedtools merge -i - -c 4 -o distinct > ${prefix}.bed

    echo "INFO: \$(wc -l < ${prefix}.bed) merged intervals from \$(wc -l < ${genes}) list entries; \$(wc -l < ${prefix}.dropped.txt) entries matched no gene" >&2
    """

    stub:
    prefix = task.ext.prefix ?: 'rescue_regions'
    """
    printf 'chr7\\t140750000\\t140760000\\tBRAF\\nchr5\\t1250000\\t1300000\\tTERT\\n' > ${prefix}.bed
    touch ${prefix}.dropped.txt
    """
}
