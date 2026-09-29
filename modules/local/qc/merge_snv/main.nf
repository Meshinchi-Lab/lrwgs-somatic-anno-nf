// modules/local/qc/merge_snv/main.nf
//
// Warning-only QC step for MERGE_SNV. Compares the final consensus VCF to the
// isec ground-truth file (0000.vcf.gz from `bcftools isec -n+2`) on:
//   1. Record count equality
//   2. CHROM:POS set equality (sorted diff)
//   3. INFO/CLAIRSTO_VAF    populated on 100% of records
//   4. INFO/CLAIRSTO_DP     populated on 100% of records
//   5. INFO/DEEPSOMATIC_VAF populated on 100% of records
//   6. INFO/DEEPSOMATIC_DP  populated on 100% of records
//
// Writes <meta.id>.merge_snv.qc.log and ALWAYS exits 0 (failures → WARN-level
// log only). Flip to a hard fail later by changing the trailing `exit 0` here
// to `exit \$status`.
//
process QC_MERGE_SNV {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(consensus_vcf), path(consensus_tbi), path(isec_0000)

    output:
    tuple val(meta), path("*.merge_snv.qc.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    LOG=${prefix}.merge_snv.qc.log
    status=0

    # Candidate-gene rescue adds single-caller records that are ABSENT from isec
    # by design (CALLER=CLAIRSTO or =DEEPSOMATIC, vs CLAIRSTO,DEEPSOMATIC for a
    # true two-caller consensus call). Comparing them against the isec ground
    # truth would warn on every rescue, so the comparison is restricted to
    # consensus-only records. The rescued count is reported separately.
    # Select consensus-only records by ABSENCE of the RESCUED tag. Do NOT test
    # INFO/CALLER=="CLAIRSTO,DEEPSOMATIC": bcftools treats a comma-separated
    # string as a value SET, so that expression matches any record containing
    # either value -- i.e. everything -- and silently filters nothing.
    bcftools view -i 'INFO/RESCUED="."' ${consensus_vcf} -Oz -o cons_only.vcf.gz
    n_rescued=\$(bcftools view -H ${consensus_vcf} | wc -l | tr -d ' ')
    n_isec=\$(bcftools view -H ${isec_0000}  | wc -l | tr -d ' ')
    n_cons=\$(bcftools view -H cons_only.vcf.gz | wc -l | tr -d ' ')
    echo "Rescued (single-caller) records: \$(( n_rescued - n_cons ))" >> \$LOG
    echo "Record count: isec=\$n_isec consensus=\$n_cons" >> \$LOG
    if [ "\$n_isec" != "\$n_cons" ]; then
        echo "WARN: record count mismatch (isec=\$n_isec consensus=\$n_cons)" >> \$LOG
        status=1
    fi

    bcftools query -f '%CHROM\\t%POS\\n' ${isec_0000}     | sort -u > sites_isec.txt
    bcftools query -f '%CHROM\\t%POS\\n' cons_only.vcf.gz | sort -u > sites_cons.txt
    if ! diff -q sites_isec.txt sites_cons.txt > /dev/null 2>&1; then
        echo "WARN: CHROM:POS set differs (see sites_*.diff)" >> \$LOG
        diff sites_isec.txt sites_cons.txt > sites_isec_vs_cons.diff || true
        status=1
    fi

    for field in CLAIRSTO_VAF CLAIRSTO_DP DEEPSOMATIC_VAF DEEPSOMATIC_DP; do
        total=\$(bcftools view -H cons_only.vcf.gz | wc -l | tr -d ' ')
        if [ "\$total" -eq 0 ]; then
            echo "INFO: \$field — 0 records to check" >> \$LOG
            continue
        fi
        # NOTE: this pattern was previously '^\\.\\?\$', which in ERE matches the
        # literal two-character string ".?" and therefore never matched anything --
        # so `grep -cv` counted every line and the check always reported 100%
        # populated regardless of the data. '^\\.\$' matches a bare "." correctly.
        populated=\$(bcftools view cons_only.vcf.gz \\
            | bcftools query -f "%INFO/\$field\\n" \\
            | grep -cvE '^\\.\$' || true)
        echo "\$field: \$populated/\$total populated" >> \$LOG
        if [ "\$populated" -ne "\$total" ]; then
            echo "WARN: \$field missing on \$((total - populated)) records" >> \$LOG
            status=1
        fi
    done

    if [ "\$status" -eq 0 ]; then
        echo "QC_MERGE_SNV: all checks passed" >> \$LOG
    else
        echo "QC_MERGE_SNV: one or more checks failed (warning-only; process exits 0)" >> \$LOG
    fi

    cat \$LOG >&2

    exit 0
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merge_snv.qc.log
    """
}
