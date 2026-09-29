// modules/local/bcftools/reheader_wakhan_fmt/main.nf
//
// Post-AnnotSV cleanup for the Wakhan CNA VCF. Two tasks:
//   1. Re-inject the Wakhan FORMAT declarations (##FORMAT=<...GT,TCN,CN1,CN2,
//      CNQ1,CNQ2,COV1,COV2>) that AnnotSV/variantconvert strips during
//      conversion. bcftools annotate --header-lines is additive — it appends
//      only, without duplicating existing keys.
//   2. Repair Wakhan's malformed BPS INFO tag. Wakhan emits:
//        ##INFO=<ID=BPS,Number=0,Type=String,Description="Breakpoints covering segment">
//      Number=0 is reserved for Type=Flag per VCF v4.2 §1.4.2, and Wakhan
//      compounds this by writing the tag inconsistently at the record level:
//      most records use proper `BPS=severus_DUP37` (key=value) but ~9% use
//      bare `BPS` (flag-only, no value). Under the malformed header cyvcf2
//      returns Python bool for the flag-only records, and vembrane's cyvcf2
//      backend then calls `.split(",")` on that bool, crashing with
//        AttributeError: 'bool' object has no attribute 'split'
//      at backend_cyvcf2.py:366 during table serialization.
//
//      The repair pipeline (auditable, bcftools-native):
//        a. Rewrite the BPS header declaration to Number=.,Type=String
//           (bcftools reheader). Records are untouched by reheader, but this
//           makes bcftools query treat BPS as a string field.
//        b. Extract BPS values via `bcftools query -f '%INFO/BPS'` into a
//           persisted TSV `${prefix}.bps_source.tsv` (emitted as an output
//           for auditability). Under the reheaded String declaration,
//           key=value records return their string; bare-flag records return
//           the literal string "1" (htslib's Flag textualization).
//        c. Normalize the "1" tokens back to "." (VCF missing-value string)
//           with a small awk pass. This is the ONLY custom text touch —
//           everything else is bcftools.
//        d. Strip the corrupt BPS field entirely from the VCF via
//           `bcftools annotate -x INFO/BPS`, then re-annotate with the
//           clean TSV via `bcftools annotate -a` + a proper
//           Number=.,Type=String header line. All records now share a
//           single, valid String encoding.
//
// KNOWN BCFTOOLS LIMITATION (upstream bug — track in a GitHub issue on
// samtools/bcftools):
//     bcftools annotate rejects records with POS=0, though 0-based
//     positions are valid per the VCF spec for telomere breakpoints
//     (chr1:0, chr21:0, etc.). Wakhan legitimately emits these.
//     Workaround: temporarily shift POS 0to1 for the target VCF + the
//     annotation source TSV before invoking bcftools annotate, then
//     shift POS 1to0 back on the annotated output. Both shifts are
//     scoped to `${id_prefix}` (currently 'wakhan') via sed on the
//     `\t<POS>\t<ID>` boundary so unrelated records are untouched.
//     Mirrors the pattern established in `bcftools/rename_ids/main.nf`
//     and `bcftools/add_caller_tag/main.nf`. Remove the shifts once
//     upstream bcftools accepts POS=0 in annotate.
//
process BCFTOOLS_REHEADER_WAKHAN_FMT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(annotsv_vcf), path(original_vcf)
    val   id_prefix                                             // e.g. 'wakhan' — scopes the POS 0↔1 sed pattern to Wakhan record IDs

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),         emit: vcf
    tuple val(meta), path("${prefix}.bps_source.tsv"), emit: bps_source     // auditable BPS extraction (CHROM POS BPS_value)
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | sed '1!d; s/^.*bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.cna.annotated"
    """
    # ── Step 1: Re-inject Wakhan FORMAT headers stripped by AnnotSV.
    bcftools view -h ${original_vcf} | grep -E '^##FORMAT' > wakhan_fmt_headers.txt

    bcftools annotate \\
        --header-lines wakhan_fmt_headers.txt \\
        -Oz \\
        -o step1.formatfix.vcf.gz \\
        ${annotsv_vcf}
    bcftools index -t step1.formatfix.vcf.gz

    # ── Step 2: Rewrite BPS declaration Number=0,Type=String to Number=.,Type=String.
    # bcftools reheader swaps header text only; records unchanged. After this
    # step, bcftools query returns BPS as strings (or "1" for flag-only records).
    bcftools view -h step1.formatfix.vcf.gz > step2_header.txt
    sed -i.bak -E 's|^##INFO=<ID=BPS,Number=0,Type=String|##INFO=<ID=BPS,Number=.,Type=String|' step2_header.txt
    bcftools reheader -h step2_header.txt -o step2.reheaded.vcf.gz step1.formatfix.vcf.gz
    bcftools index -t step2.reheaded.vcf.gz

    # ── Step 3: Extract BPS values via bcftools query to auditable TSV.
    # Under the reheaded String declaration, key=value records return their
    # real severus IDs; bare-flag records return the literal string "1".
    # A one-liner awk normalizes "1" to "." (VCF missing-value string).
    bcftools query -f '%CHROM\\t%POS\\t%INFO/BPS\\n' step2.reheaded.vcf.gz \\
        | awk -F '\\t' 'BEGIN{OFS="\\t"} { if (\$3 == "1") \$3 = "."; print }' \\
        > ${prefix}.bps_source.tsv

    # ── Step 4: Strip corrupt BPS field entirely — clears mixed encoding.
    bcftools annotate -x INFO/BPS -Oz -o step4.stripped.vcf.gz step2.reheaded.vcf.gz
    bcftools index -t step4.stripped.vcf.gz

    # ── Step 5: Re-inject BPS via bcftools annotate with proper Number=.,Type=String.
    cat > bps_hdr.txt <<EOF
##INFO=<ID=BPS,Number=.,Type=String,Description="Breakpoints covering segment (severus IDs from Wakhan)">
EOF

    # WORKAROUND for bcftools annotate not accepting POS=0 (see module header comment):
    # shift POS 0to1 in both the annotation source TSV and the target VCF for
    # the annotate call, then shift 1to0 back on the output. Scoped to
    # ${id_prefix}-prefixed IDs so non-Wakhan records are unaffected.

    # 5a. Shift POS 0to1 in the TSV, then bgzip + tabix.
    awk -F '\\t' 'BEGIN{OFS="\\t"} { if (\$2 == "0") \$2 = "1"; print }' \\
        ${prefix}.bps_source.tsv \\
        | bgzip -c > ${prefix}.bps_source.1pos.tsv.gz
    tabix -s1 -b2 -e2 ${prefix}.bps_source.1pos.tsv.gz

    # 5b. Shift POS 0to1 in the target VCF (scoped to ${id_prefix} IDs).
    bcftools view step4.stripped.vcf.gz \\
        | sed -E "s|\\t0\\t${id_prefix}|\\t1\\t${id_prefix}|" > step4.1pos.vcf

    # 5c. Annotate; revert POS 1to0 for ${id_prefix} IDs; recompress + index.
    bcftools annotate \\
        -a ${prefix}.bps_source.1pos.tsv.gz \\
        -h bps_hdr.txt \\
        -c CHROM,POS,INFO/BPS \\
        step4.1pos.vcf \\
        | sed -E "s|\\t1\\t${id_prefix}|\\t0\\t${id_prefix}|" \\
        | bcftools view -Oz -o ${prefix}.vcf.gz
    bcftools index -t ${prefix}.vcf.gz

    # ── Cleanup intermediates (retain ${prefix}.bps_source.tsv for auditability).
    rm -f step1.formatfix.vcf.gz step1.formatfix.vcf.gz.tbi \\
          step2.reheaded.vcf.gz step2.reheaded.vcf.gz.tbi \\
          step2_header.txt step2_header.txt.bak \\
          step4.stripped.vcf.gz step4.stripped.vcf.gz.tbi step4.1pos.vcf \\
          wakhan_fmt_headers.txt bps_hdr.txt \\
          ${prefix}.bps_source.1pos.tsv.gz ${prefix}.bps_source.1pos.tsv.gz.tbi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.cna.annotated"
    """
    echo '' | gzip > ${prefix}.vcf.gz
    printf 'CHROM\\tPOS\\tBPS\\n' > ${prefix}.bps_source.tsv
    """
}
