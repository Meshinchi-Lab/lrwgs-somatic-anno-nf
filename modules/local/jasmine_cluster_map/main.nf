// Turn Jasmine's merged VCF into a long, per-(sample x cluster) recurrence TSV.
//
// Replaces the reference repo's transfer_annotations.py. Annotations are never copied INTO
// the merged VCF; instead the cluster map is joined back out onto tables that already hold
// them. Every VCF field is extracted natively by `bcftools query` -- the awk below is
// relational TSV work, not VCF parsing.
//
// THE JOIN IS OUTER, and that is the point. Jasmine clusters the PRE-FILTER VCFs, so a
// cluster can contain records absent from the tiered report TSV -- sub-threshold calls in
// samples that filtered the SV out. Those are the rescue candidates the whole subworkflow
// exists to surface; an inner join would silently drop them. Each output row is labelled
// TIERED or RESCUED, and rescued rows take their annotations from the record table.
//
// Two counts express the outcome:
//   n_tiered = samples where this SV was a Tier 1/2/3 call
//   n_any    = samples with ANY evidence in the cluster (Jasmine's INFO/SUPP)
// n_tiered=1, n_any=5 -> likely false negative (over-filtered elsewhere)
// n_tiered=1, n_any=1 -> likely false positive (no corroboration)
process JASMINE_CLUSTER_MAP {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  report_tsv      // cohort_sv_pass2.tsv -- the TIERED set (driver_tier lives here)
    path  records_tsv     // cohort-wide bcftools query of ALL in-region records
    val   min_recurrence
    val   max_dist        // only used to flag suspiciously wide (chained) clusters

    output:
    tuple val(meta), path("${prefix}.tsv"),      emit: tsv
    tuple val(meta), path("${prefix}.map.tsv"),  emit: map
    tuple val("${task.process}"), val('bcftools'),
          eval("bcftools --version | head -1 | sed 's/bcftools //'"),
          topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.cohort_sv_recurrent"
    """
    # 1. Cluster map, native bcftools. One row per (cluster, source_id).
    #    IDLIST is comma-separated; SUPP_VEC is a binary string, one char per input sample.
    #    INFO/SUPP is emitted natively (Number=1,Type=Integer) and equals the set-bit count
    #    of SUPP_VEC, so no bit-counting is needed.
    bcftools query -f '%ID\\t%INFO/IDLIST\\t%INFO/SUPP_VEC\\t%INFO/SUPP\\n' ${vcf} \\
      | awk -F'\\t' 'BEGIN{OFS="\\t"; print "cluster_id","source_id","supp_vec","n_any"}
          {
              n = split(\$2, ids, ",")
              for (i=1; i<=n; i++) print \$1, ids[i], \$3, \$4
          }' > ${prefix}.map.tsv

    # 2. Three-way OUTER join, uniform schema for TIERED and RESCUED rows.
    #    awk is used because coreutils `join` cannot express a 3-way join with a computed
    #    per-group aggregate (n_tiered, cluster_span).
    awk -F'\\t' -v OFS='\\t' -v minrec=${min_recurrence} \\
        -v mapf='${prefix}.map.tsv' -v repf='${report_tsv}' '
        # pass 1 -- cluster membership
        FILENAME==mapf { if (FNR>1) { cl[\$2]=\$1; nany[\$1]=\$4 } ; next }

        # pass 2 -- the TIERED set, keyed SAMPLE_ID to match the uniquified source_id
        FILENAME==repf {
            if (FNR==1) { for (i=1;i<=NF;i++) c[\$i]=i ; next }
            tier[ \$(c["SAMPLE"]) "_" \$(c["ID"]) ] = \$(c["driver_tier"])
            next
        }

        # pass 3 -- every in-region record (bcftools query output, no header)
        {
            src = \$2
            if (!(src in cl)) { orphan++; next }   # expect 0: Jasmine clusters singletons too
            k = cl[src]
            line[src] = \$0 ; clust[src] = k ; order[++nrec] = src
            if (src in tier && !((k SUBSEP \$1) in seen)) { seen[k SUBSEP \$1]=1 ; nt[k]++ }
            # Observed cluster width. Jasmine merges by SINGLE LINKAGE, so a cluster can be
            # far wider than max_dist when intermediates bridge the extremes. This does not
            # prevent chaining -- it makes it visible.
            pos = \$4 + 0
            if (!(k in mn) || pos < mn[k]) mn[k] = pos
            if (!(k in mx) || pos > mx[k]) mx[k] = pos
        }

        END {
            print "cluster_id","source_id","sample","chrom","pos","end","svtype","svlen",
                  "gene_name","gnomad_sv_af","clinvar_clnsig","caller",
                  "record_class","driver_tier","n_tiered","n_any","cluster_span","cohort_recurrent"
            for (i=1; i<=nrec; i++) {
                src = order[i] ; k = clust[src]
                split(line[src], f, "\\t")
                rc = (src in tier) ? "TIERED"   : "RESCUED"
                dt = (src in tier) ? tier[src]  : "NA"
                nT = (k in nt)     ? nt[k]      : 0
                print k, src, f[1], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11],
                      rc, dt, nT, nany[k], (mx[k] - mn[k]), (nT >= minrec ? "TRUE" : "FALSE")
            }
            if (orphan > 0)
                printf("WARN: %d in-region records matched no Jasmine cluster\\n", orphan) > "/dev/stderr"
        }' ${prefix}.map.tsv ${report_tsv} ${records_tsv} > ${prefix}.tsv

    # 3. QC -- the dual counts ARE the result, so report their distribution.
    awk -F'\\t' -v OFS='\\t' '
        NR==1 { for (i=1;i<=NF;i++) c[\$i]=i ; next }
        {
            total++
            if (\$(c["record_class"])=="TIERED")  tiered++ ; else rescued++
            if (\$(c["cohort_recurrent"])=="TRUE") rec++
            nT=\$(c["n_tiered"]) ; nA=\$(c["n_any"])
            if (nT==1 && nA>=2) fn++          # tiered once, evidence elsewhere -> false negative
            if (nT==1 && nA==1) fp++          # tiered once, no corroboration   -> false positive
            if (\$(c["cluster_span"]) + 0 > 2 * ${max_dist}) wide++
        }
        END {
            printf("INFO: %d rows (%d TIERED, %d RESCUED); %d recurrent (n_tiered>=%s)\\n",
                   total, tiered, rescued, rec, "${min_recurrence}") > "/dev/stderr"
            printf("INFO: %d rescue candidates (n_tiered=1, n_any>=2); %d private/FP candidates (n_tiered=1, n_any=1)\\n",
                   fn, fp) > "/dev/stderr"
            if (wide > 0)
                printf("WARN: %d rows sit in clusters wider than 2x max_dist (%d bp) -- single-linkage chaining; inspect cluster_span\\n",
                       wide, 2 * ${max_dist}) > "/dev/stderr"
        }' ${prefix}.tsv

    # 3b. Guard: Jasmine can leave INFO/END holding a MATE coordinate on another chromosome,
    #     giving END < POS (Jasmine issue #26). htslib only WARNS and silently drops the tag
    #     from the index, so this must be checked explicitly. Measured as 0 for this
    #     pipeline's inputs; the reference implementation
    #     (jasmine_sv_annotations/scripts/clean_bnds.py, Hiatt 2024) exists solely to repair
    #     it. If this fires, repair it natively with
    #     `bcftools annotate -a <tsv> -c CHROM,POS,INFO/END` setting END=POS for BND rows.
    badend=\$(bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\n' ${vcf} \\
              | awk -F'\\t' '\$3!="." && \$3+0 < \$2+0' | wc -l | tr -d ' ')
    if [ "\$badend" -gt 0 ]; then
        echo "WARN: \$badend merged records have INFO/END < POS (Jasmine issue #26); htslib will ignore the END tag when indexing" >&2
    fi

    # 4. Every TIERED key must have landed in a cluster; a miss means the source_id prefix
    #    no longer matches SAMPLE and the join key is broken.
    missing=\$(awk -F'\\t' '
        FNR==NR { if (FNR>1) have[\$2]=1 ; next }
        FNR==1  { for (i=1;i<=NF;i++) c[\$i]=i ; next }
        !( (\$(c["SAMPLE"]) "_" \$(c["ID"])) in have ) { m++ }
        END { print m+0 }' ${prefix}.map.tsv ${report_tsv})
    if [ "\$missing" -gt 0 ]; then
        echo "WARN: \$missing tiered rows had no cluster -- check COHORT_SV_SET_ID's id prefix matches SAMPLE" >&2
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.cohort_sv_recurrent"
    """
    printf 'cluster_id\\tsource_id\\tsupp_vec\\tn_any\\n' > ${prefix}.map.tsv
    printf 'cluster_id\\tsource_id\\tsample\\tchrom\\tpos\\tend\\tsvtype\\tsvlen\\tgene_name\\tgnomad_sv_af\\tclinvar_clnsig\\tcaller\\trecord_class\\tdriver_tier\\tn_tiered\\tn_any\\tcluster_span\\tcohort_recurrent\\n' > ${prefix}.tsv
    """
}
