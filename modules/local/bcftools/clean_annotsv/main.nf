// modules/local/bcftools/clean_annotsv/main.nf
//
// Normalize AnnotSV's compound-encoded ACMG_class + AnnotSV_ranking_score to
// clean scalars BEFORE FILTER_SV / FILTER_CNA so downstream vembrane, KnotAnnotSV,
// and the R report see them as raw Number=1 typed fields instead of Python-tuple
// wrapped Number=. strings.
//
// AnnotSV emits ACMG_class as `<split_mode_value>,full=<full_mode_value>` on
// records that carry Annotation_mode=full,split (see raw output — e.g. `3,full=3`,
// `.,full=NA`, `1,full=1`). The clinical interpretation uses the `full=X` value
// (the full-region assessment). This module extracts that value into a bare
// integer + rewrites the header from Number=.,Type=String to Number=1,Type=Integer.
//
// AnnotSV_ranking_score is already a clean integer in the raw VCF; only its
// header type needs correcting (Number=1,Type=Float).
//
// Mechanism (three-step bcftools chain — no custom scripts outside the module):
//   1. `bcftools query` extracts CHROM/POS/REF/ALT + the two INFO fields; awk
//      re-encodes ACMG_class to just the full-mode integer (NA -> `.`).
//   2. `bcftools annotate --remove` strips the old Number=.,Type=String headers
//      and values.
//   3. `bcftools annotate --header-lines --annotations` re-adds the two fields
//      with the corrected Number=1 typed header lines from the cleaned TSV.
//
// Inputs:
//   tuple val(meta), path(vcf)  — output of BCFTOOLS_REHEADER_*_FMT (no tbi)
//
// Outputs:
//   tuple val(meta), path("*.acmg_clean.vcf.gz"), emit: vcf
//

process BCFTOOLS_CLEAN_ANNOTSV {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.acmg_clean.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.acmg_clean"
    """
    # Step 1: extract clean scalars into a bgzipped + tabix-indexed TSV.
    # ACMG_class raw form is `<split>,full=<full>` or `.,full=NA` — pull out
    # the full=X value; NA / empty / '.' -> missing (`.`).
    # AnnotSV_ranking_score is already a clean integer; coerce empty -> `.`.
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/ACMG_class\\t%INFO/AnnotSV_ranking_score\\n' \\
        ${vcf} \\
        | awk -F'\\t' -v OFS='\\t' '{
            acmg = "."
            if (\$5 != "" && \$5 != ".") {
                n = split(\$5, parts, ",")
                for (i = 1; i <= n; i++) {
                    if (parts[i] ~ /^full=/) {
                        val = parts[i]
                        sub(/^full=/, "", val)
                        if (val != "NA" && val != "" && val != ".") acmg = val
                        break
                    }
                }
                # Fallback: no full=X prefix seen, but value is a bare integer
                if (acmg == "." && \$5 ~ /^[0-9]+\$/) acmg = \$5
            }
            score = (\$6 == "" || \$6 == ".") ? "." : \$6
            print \$1, \$2, \$3, \$4, acmg, score
        }' \\
        | bgzip -c > clean.tsv.gz
    tabix -s1 -b2 -e2 clean.tsv.gz

    # Step 2: drop the old Number=.,Type=String headers + values so re-annotate
    # in step 3 uses OUR new header definitions instead of re-using the old ones.
    # bcftools annotate --header-lines is additive-only — it will silently skip
    # duplicate INFO IDs already declared in the target, so the --remove pass is
    # mandatory to actually replace the type.
    bcftools annotate \\
        --remove 'INFO/ACMG_class,INFO/AnnotSV_ranking_score' \\
        ${vcf} \\
        -Oz -o stripped.vcf.gz

    # Step 3: re-add both fields with clean scalar header definitions.
    cat > new_hdr.txt <<'EOF'
##INFO=<ID=ACMG_class,Number=1,Type=Integer,Description="AnnotSV full-mode ACMG classification (1-5); missing when NA">
##INFO=<ID=AnnotSV_ranking_score,Number=1,Type=Float,Description="AnnotSV ranking score">
EOF

    bcftools annotate \\
        --header-lines new_hdr.txt \\
        -a clean.tsv.gz \\
        -c CHROM,POS,REF,ALT,INFO/ACMG_class,INFO/AnnotSV_ranking_score \\
        stripped.vcf.gz \\
        -Oz -o ${prefix}.vcf.gz

    # Cleanup intermediates (Nextflow catches versions.yml only)
    rm -f clean.tsv.gz clean.tsv.gz.tbi stripped.vcf.gz new_hdr.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.acmg_clean"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    """
}
