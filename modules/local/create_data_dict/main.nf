// modules/local/create_data_dict/main.nf
//
// Parse a VCF's ##INFO / ##FORMAT / ##FILTER header lines into a TSV
// "data dictionary" (one row per declared field). Consumed downstream by
// SV_REPORT_INDEX to render definition tables above each DT table in the
// cohort report, so reviewers see the source callers' own descriptions
// rather than a downstream author's paraphrase.
//
// Columns emitted (tab-separated, one header row):
//   field_type      INFO | FORMAT | FILTER
//   field_name      ID= attribute
//   number          Number= attribute (blank for FILTER lines — n/a)
//   type            Type= attribute (blank for FILTER lines — n/a)
//   description     Description= attribute, quotes stripped, tabs replaced
//                   with single spaces so the TSV stays parseable
//   source          Source= attribute (blank when absent)
//
// Implementation: `bcftools view -h` yields the full header; awk parses each
// line by explicitly locating each Key=<value> attribute via POSIX-safe
// regex primitives (`match()` + RSTART/RLENGTH), handling quoted Description
// values (which may embed commas) as an atomic unit. No external Python
// dependency needed — runs in the same bcftools_htslib container the rest
// of the CNA branch uses.
//
// Runs once per sample, once per variant type (SV, CNA, SNV). The three
// FILTER_{SV,CNA,SNV} subworkflows call the process aliased differently.
//
process CREATE_DATA_DICT {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.data_dict.tsv"), emit: dict
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | sed '1!d; s/^.*bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # `bcftools view -h` handles both .vcf.gz and .vcf transparently.
    # The awk program below is POSIX-compliant (uses match()+RSTART/RLENGTH,
    # no gawk-only third-arg to match()); works in mawk/busybox awk too.
    bcftools view -h ${vcf} \\
        | awk 'BEGIN {
            OFS = "\\t"
            print "field_type", "field_name", "number", "type", "description", "source"
        }
        /^##(INFO|FORMAT|FILTER)=</ {
            line = \$0

            # field_type: characters between ^## and =<
            ftype = line
            sub(/=<.*\$/, "", ftype)
            sub(/^##/,   "", ftype)

            # content: characters inside <...>
            content = line
            sub(/^##[A-Z]+=</, "", content)
            sub(/>\$/,          "", content)

            id  = ""; num = ""; typ = ""; desc = ""; src = ""

            # ID= up to next comma or end
            work = content
            if (match(work, /ID=[^,>]+/)) {
                id = substr(work, RSTART + 3, RLENGTH - 3)
            }

            # Number= (only meaningful for INFO/FORMAT; blank for FILTER)
            work = content
            if (match(work, /Number=[^,>]+/)) {
                num = substr(work, RSTART + 7, RLENGTH - 7)
            }

            # Type= (only meaningful for INFO/FORMAT; blank for FILTER)
            work = content
            if (match(work, /Type=[^,>]+/)) {
                typ = substr(work, RSTART + 5, RLENGTH - 5)
            }

            # Description= is quoted per VCF spec and may embed commas.
            # Match the full quoted value first; fall back to unquoted for
            # any non-conforming producers.
            work = content
            if (match(work, /Description="[^"]*"/)) {
                # +13 skips  Description="  ; -14 drops opening + closing "
                desc = substr(work, RSTART + 13, RLENGTH - 14)
            } else if (match(work, /Description=[^,>]+/)) {
                desc = substr(work, RSTART + 12, RLENGTH - 12)
            }

            # Source= (optional; often absent for VCF v4.2 core fields)
            work = content
            if (match(work, /Source="[^"]*"/)) {
                src = substr(work, RSTART + 8, RLENGTH - 9)
            } else if (match(work, /Source=[^,>]+/)) {
                src = substr(work, RSTART + 7, RLENGTH - 7)
            }

            # Guard the TSV: strip any literal tabs or newlines that a producer
            # smuggled into a Description value — either would break the TSV
            # row structure downstream.
            gsub(/\\t/, " ", desc); gsub(/\\r/, " ", desc); gsub(/\\n/, " ", desc)
            gsub(/\\t/, " ", src);  gsub(/\\r/, " ", src);  gsub(/\\n/, " ", src)

            print ftype, id, num, typ, desc, src
        }' > ${prefix}.data_dict.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'field_type\\tfield_name\\tnumber\\ttype\\tdescription\\tsource\\n' \\
        > ${prefix}.data_dict.tsv
    """
}
