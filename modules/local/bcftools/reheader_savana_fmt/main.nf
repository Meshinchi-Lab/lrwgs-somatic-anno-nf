// Post-AnnotSV header repair for the SV VCF, in one bcftools annotate call:
//   1. Re-inject the Severus/SAVANA FORMAT declarations (VAF, hVAF, DR, DV). AnnotSV does not
//      carry the input's ##FORMAT lines over, so these fields are present in records but
//      undefined in the header, and bcftools / pysam / vembrane reject or mistype them.
//      bcftools annotate --header-lines is additive: it appends only.
//   2. Put contig names back on the reference FASTA's. AnnotSV strips "chr" from CHROM, declares
//      its own Ensembl-style contigs (1..22, X, Y, MT and ~170 KI/GL accession contigs; it writes records only
//      on 1..22, X, Y and M/MT) and writes a chrM call as "M" under that
//      "MT" header, so bcftools stopped on any sample with a chrM call ("Contig 'M' is not
//      defined in the header") and FASTA lookups downstream failed (Jasmine --dup_to_ins asked
//      samtools faidx for "11" from a FASTA holding "chr11"). The --rename-chrs map and the final
//      ##contig lines come from the FASTA's .fai, so a UCSC-style reference gives chr1..chrM and
//      an Ensembl-style one 1..MT. The ##contig=<ID=M> header line lets bcftools read the "M"
//      record before it is renamed. Assumes params.fasta is the calling reference, so BND mates
//      in ALT (copied unchanged by AnnotSV) already carry its names. A CHROM with no FASTA contig
//      fails the task with the names involved.
// BCFTOOLS_REHEADER_WAKHAN_FMT applies the same rename to the CNA VCF.
process BCFTOOLS_REHEADER_SAVANA_FMT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fai)     // samtools faidx index of params.fasta: the contig-name authority

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.ensemble"

    """
    # hVAF is declared Number=3,Type=Float (HP1, HP2, unphased) — SEVERUS emits three comma-separated values per record.
    # ##contig=<ID=M> declares AnnotSV's mitochondrial CHROM so bcftools can read it before the rename.
    cat > savana_fmt_headers.txt << 'HDRS'
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Tumour variant allele frequency">
##FORMAT=<ID=hVAF,Number=3,Type=Float,Description="Haplotype-aware tumour variant allele frequency (HP1,HP2,unphased ONT haplotag)">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Tumour reference-supporting reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Tumour variant-supporting reads">
##contig=<ID=M>
HDRS

    # Rename map from the FASTA index: each "chr" contig from its bare name (1 -> chr1, M -> chrM).
    # AnnotSV writes the mitochondrial CHROM as M, so a FASTA whose mitochondrial contig is MT or
    # chrMT also needs M -> that name. Names AnnotSV leaves unchanged are not in the map.
    cut -f1 ${fai} | sed -nE 's/^chr(.+)\$/\\1\\tchr\\1/p' > contig_map.tsv
    mito=\$(cut -f1 ${fai} | grep -m1 -xE '(chr)?(M|MT)' || true)
    case "\$mito" in MT|chrMT) printf 'M\\t%s\\n' "\$mito" >> contig_map.tsv ;; esac

    bcftools annotate \\
        --header-lines savana_fmt_headers.txt \\
        --rename-chrs contig_map.tsv \\
        -Oz \\
        -o renamed.vcf.gz \\
        ${vcf}

    # Unmapped check: every CHROM must now be a FASTA contig, or its records could not be looked
    # up in the FASTA downstream. bcftools does not stop on this by itself.
    unmapped=\$(bcftools query -f '%CHROM\\n' renamed.vcf.gz | LC_ALL=C sort -u \\
        | { grep -vxF -f <(cut -f1 ${fai}) || [ \$? -eq 1 ]; } | paste -sd, -)
    if [ -n "\$unmapped" ]; then
        echo "ERROR: ${vcf} has CHROM values with no contig in the reference FASTA index ${fai}: \$unmapped" >&2
        exit 1
    fi

    # Header contigs = the FASTA's names, lengths and order. reheader --fai keeps the order of the
    # contig lines already present (AnnotSV's lexical 1, 10, 11, ...), which every downstream
    # bcftools sort follows, so drop them first and let it add the .fai's in .fai order.
    bcftools view -h renamed.vcf.gz | grep -v '^##contig=' > header_no_contigs.txt
    bcftools reheader -h header_no_contigs.txt -o no_contigs.vcf.gz renamed.vcf.gz
    bcftools reheader --fai ${fai} -o ${prefix}.vcf.gz no_contigs.vcf.gz
    rm -f renamed.vcf.gz no_contigs.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.ensemble"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    """
}
