// modules/local/setupclinvar/main.nf
//
//
// Inputs:
//   path(clinvar_vcf)  — raw ClinVar VCF (bgzipped). Index (.tbi) must be at
//                         <clinvar_vcf>.tbi (Nextflow stages siblings via
//                         path(); we pass them together).
//
// Outputs:
//   path("clinvar.norm.vcf.gz")     — split-multiallelic, sorted, indexed
//   path("clinvar.norm.vcf.gz.tbi")  — tabix index
//

process SETUPCLINVAR {
    tag "Clinical VCF GRCh38"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    path clinvar_vcf
    path clinvar_tbi   // staged alongside; not used by name but required for tabix to find it
    path rename_annots

    output:
    path "*.norm.vcf.gz",     emit: vcf
    path "*.norm.vcf.gz.tbi", emit: tbi
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: 'clinvar'
    rename_annots = rename_annots ? "--rename-annots ${rename_annots}" : ''
    """
    # Introduced in 2024 somatic-specific tags (SCI / ONC). 
    have_sci=\$(bcftools view -h ${clinvar_vcf} | grep -c "^##INFO=<ID=SCI," || true)
    have_onc=\$(bcftools view -h ${clinvar_vcf} | grep -c "^##INFO=<ID=ONC," || true)
    if [ "\$have_sci" -eq 0 ] || [ "\$have_onc" -eq 0 ]; then
        echo "WARN: SETUPCLINVAR — ClinVar VCF predates 2024 restructure (SCI=\$have_sci ONC=\$have_onc)." >&2
    else
        echo "INFO: SETUPCLINVAR — 2024 SCI/ONC tags present in header." >&2
    fi

    # Chromosome-naming reconciliation: raw ClinVar from NCBI uses unprefixed
    # contigs (1..22, X, Y, MT) but the pipeline FASTA (GCA_000001405.15
    # GRCh38_no_alt_analysis_set) uses UCSC-style chr-prefixed contigs
    # (chr1..chr22, chrX, chrY, chrM). 
    #
    # Build the rename map from ClinVar's header — UCSC convention: MT -> chrM.
    bcftools view -h ${clinvar_vcf} \\
        | grep '^##contig' \\
        | sed -E 's/^##contig=<ID=([^,>]+).*/\\1/' \\
        | awk -v OFS='\\t' '
            /^chr/                          { next }
            /^([1-9]|1[0-9]|2[0-2])\$/      { print \$0, "chr" \$0; next }
            \$0 == "X" || \$0 == "Y"        { print \$0, "chr" \$0; next }
            \$0 == "MT"                     { print \$0, "chrM";     next }
        ' > chr_rename.tsv

    if [ -s chr_rename.tsv ]; then
        echo "INFO: SETUPCLINVAR — renaming \$(wc -l < chr_rename.tsv) unprefixed ClinVar contigs to chr-prefixed UCSC style." >&2
        bcftools annotate \\
            ${rename_annots} \\
            --rename-chrs chr_rename.tsv \\
            ${clinvar_vcf} -Oz \\
            -o ${prefix}.renamed.vcf.gz
        CLINVAR_IN=${prefix}.renamed.vcf.gz
    else
        echo "INFO: SETUPCLINVAR — ClinVar contigs already chr-prefixed; skipping rename." >&2
        CLINVAR_IN=${clinvar_vcf}
    fi

    # Normalize: split multi-allelics (-m -both), keep ClinVar's existing
    # left-alignment as-is (no -f reference needed — ClinVar is pre-aligned).
    bcftools norm -m -both "\$CLINVAR_IN" -Oz -o ${prefix}.norm.unsorted.vcf.gz

    # Sort + index the normalized output.
    bcftools sort ${prefix}.norm.unsorted.vcf.gz -Oz -o ${prefix}.norm.vcf.gz
    bcftools index -t ${prefix}.norm.vcf.gz

    rm -f *.norm.unsorted.vcf.gz *.renamed.vcf.gz chr_rename.tsv
    """

    stub:
    prefix = task.ext.prefix ?: 'clinvar'
    """
    touch ${prefix}.norm.vcf.gz
    touch ${prefix}.norm.vcf.gz.tbi
    """
}
