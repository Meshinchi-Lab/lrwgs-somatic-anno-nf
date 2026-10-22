// modules/local/sv_reference_filter_check/main.nf
//
// FILTER-column census of a VEP --custom reference VCF (gnomAD-SV sites,
// ClinVar/dbVar nstd102 SV) and, by policy, a PASS-only subset for VEP.
//
// Why: `--custom ... type=overlap` reports EVERY reference record that spans a
// query position, whatever its FILTER column says. gnomAD-SV v4.1 sites carries
// 955,369 non-PASS records out of 2,154,486 (UNRESOLVED, LOWQUAL_WHAM_SR_DEL,
// OUTLIER_SAMPLE_ENRICHED, HIGH_NCR, FAIL_MANUAL_REVIEW, ...; measured
// 2026-10-22) -- allele frequencies gnomAD itself does not stand behind -- and
// max(gnomAD_SV_AF) is a FILTER_SV gate. dbVar's nstd102 has no FILTER
// annotation at all ('.'); that is "unfiltered", not "filtered", and is
// accepted as-is.
//
// Policy (params.sv_reference_filter_policy):
//   filter  (default) non-PASS / REMOVE-labelled records are dropped
//                     (`bcftools view -f PASS,.`) and the re-indexed subset
//                     replaces the reference for VEP. A clean file is passed
//                     through untouched (symlink, no copy).
//   fail              abort the run when any such record exists.
//   warn              report only; VEP annotates against the file as supplied.
// Every policy writes <stem>.filter_check.tsv (one-line verdict) and
// <stem>.filter_counts.tsv (FILTER value census), published under
// reference/sv_filter_check/.
//
// "REMOVE": none of the references used here carries such a label (checked
// 2026-10-22); the FILTER + INFO scan exists so a future reference that marks
// records for removal trips the same check.
//
// The checked file keeps its original basename under checked/ so the
// `--custom file=<basename>` argument built in modules.config is unchanged.
process SV_REFERENCE_FILTER_CHECK {
    tag "${vcf.name}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple path(vcf), path(tbi)
    val   policy

    output:
    tuple path("checked/${vcf.name}"), path("checked/${vcf.name}.tbi"), emit: vcf_tbi
    path "*.filter_check.tsv",  emit: summary
    path "*.filter_counts.tsv", emit: counts
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def name = vcf.name
    def stem = name.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
    if ( !(policy in ['filter', 'fail', 'warn']) ) {
        error "SV_REFERENCE_FILTER_CHECK: sv_reference_filter_policy must be one of filter | fail | warn, got '${policy}'"
    }
    """
    mkdir -p checked

    # 1. FILTER census. "." = no filter applied (dbVar); PASS = passed; anything
    #    else is a filter flag (several are ';'-joined on one record).
    bcftools query -f '%FILTER\\n' ${vcf} \\
      | sort | uniq -c \\
      | awk 'BEGIN{OFS="\\t"; print "filter","n"} {print \$2, \$1}' \\
      > ${stem}.filter_counts.tsv

    total=\$(awk -F'\\t' 'NR>1 {s+=\$2} END{print s+0}' ${stem}.filter_counts.tsv)
    pass=\$(awk -F'\\t' 'NR>1 && \$1=="PASS" {s+=\$2} END{print s+0}' ${stem}.filter_counts.tsv)
    unfiltered=\$(awk -F'\\t' 'NR>1 && \$1=="." {s+=\$2} END{print s+0}' ${stem}.filter_counts.tsv)
    filtered=\$(( total - pass - unfiltered ))

    # 2. REMOVE-label scan over FILTER + INFO (whole word, case-insensitive).
    remove=\$(bcftools view -H ${vcf} | cut -f7,8 | grep -c -i -w 'REMOVE' || true)

    status=OK
    if [ "\$filtered" -gt 0 ] || [ "\$remove" -gt 0 ]; then status=NON_PASS_PRESENT; fi

    # 3. Policy.
    action=PASSTHROUGH
    kept=\$total
    case "${policy}" in
      filter)
        if [ "\$status" != OK ]; then
          bcftools view -f PASS,. --threads ${task.cpus} -Oz -o checked/${name} ${vcf}
          bcftools index -t checked/${name}
          kept=\$(bcftools view -H checked/${name} | wc -l | tr -d ' ')
          action=FILTERED_TO_PASS
        fi ;;
      fail)
        if [ "\$status" != OK ]; then
          echo "ERROR: SV_REFERENCE_FILTER_CHECK -- ${name}: \$filtered of \$total records carry a non-PASS FILTER and \$remove are labelled REMOVE. Supply a PASS-only reference (bcftools view -f PASS,. -Oz -o <out> <in> && bcftools index -t <out>) or set sv_reference_filter_policy to 'filter' or 'warn'." >&2
          exit 1
        fi ;;
      warn)
        if [ "\$status" != OK ]; then
          echo "WARN: SV_REFERENCE_FILTER_CHECK -- ${name}: \$filtered of \$total records carry a non-PASS FILTER and \$remove are labelled REMOVE; VEP annotates against the file AS SUPPLIED (sv_reference_filter_policy=warn)." >&2
        fi ;;
    esac
    if [ "\$action" = PASSTHROUGH ]; then
      ln -s "\$(readlink -f ${vcf})" checked/${name}
      ln -s "\$(readlink -f ${tbi})" checked/${name}.tbi
    fi

    # 4. One-line verdict.
    printf 'file\\ttotal\\tpass\\tunfiltered\\tfiltered\\tremove_labelled\\tstatus\\tpolicy\\taction\\trecords_for_vep\\n' \\
      > ${stem}.filter_check.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
      "${name}" "\$total" "\$pass" "\$unfiltered" "\$filtered" "\$remove" "\$status" "${policy}" "\$action" "\$kept" \\
      >> ${stem}.filter_check.tsv
    echo "INFO: SV_REFERENCE_FILTER_CHECK -- ${name}: total=\$total PASS=\$pass unfiltered=\$unfiltered filtered=\$filtered remove=\$remove status=\$status policy=${policy} action=\$action records_for_vep=\$kept" >&2
    """

    stub:
    def name = vcf.name
    def stem = name.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
    """
    mkdir -p checked
    ln -s "\$(readlink -f ${vcf})" checked/${name}
    ln -s "\$(readlink -f ${tbi})" checked/${name}.tbi
    printf 'filter\\tn\\nPASS\\t0\\n' > ${stem}.filter_counts.tsv
    printf 'file\\ttotal\\tpass\\tunfiltered\\tfiltered\\tremove_labelled\\tstatus\\tpolicy\\taction\\trecords_for_vep\\n%s\\t0\\t0\\t0\\t0\\t0\\tOK\\t%s\\tPASSTHROUGH\\t0\\n' "${name}" "${policy}" > ${stem}.filter_check.tsv
    """
}
