// Turn a cleaned report TSV into a BED of regions to restrict cohort merging to.
//
// Variant-type agnostic: parts 2 (CNA) and 3 (SNV) of the cohort merge reuse this
// unchanged by passing different column names.
//
// PADDING IS MANDATORY, not cosmetic. The whole purpose of the cohort merge is to
// cluster the SAME event called at DIFFERENT coordinates in different samples. A BED
// built from one sample's exact breakpoints would exclude another sample's copy of that
// event before the clusterer ever sees it -- silently defeating the feature. The pad must
// be >= the clusterer's distance tolerance (params.cohort_sv_max_dist).
//
// CONTIG NAMING: the cleaned TSV and the annotated VCFs both use UNPREFIXED contigs
// ("14", not "chr14"), so no renaming is applied here. Verified 2026-09-17. If the
// report layer ever starts emitting chr-prefixed CHROM, this module must be revisited.
//
// BND ROWS: a breakend is a POINT, and its interval must be POS +/- pad -- never POS..END.
//
// For INTERchromosomal BNDs the TSV carries END == POS and SVLEN == 0, so an END-based
// span would be harmless. But for INTRAchromosomal BNDs END is the MATE coordinate on the
// same contig, and using it produces a span instead of a breakpoint. Measured on the
// 14-row test set (2026-09-17): one such row,
//     9  POS=128694156  END=131151786  SVTYPE=BND  tier=1
// yielded a 2,458,031 bp interval -- 67% of the entire region set from a single record.
// At cohort scale this drags unrelated SVs into every such locus and inflates the
// clusterer's input. Hence the explicit SVTYPE=="BND" guard below.
//
// The mate is NOT parsed out of ALT (e.g. "N[chr5:171319244[", where the mate IS
// chr-prefixed while CHROM is not): each breakend is emitted by the caller as its OWN VCF
// record with its own CHROM/POS, so it already has its own TSV row and its own interval.
// Parsing ALT would duplicate intervals and require stripping the inconsistent prefix.
process COLLECT_COORDINATES {
    tag "${tsv.simpleName}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(tsv)
    val  pad
    val  tiers          // comma-separated driver_tier values to keep, e.g. '1,2,3'

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'),
          eval("bedtools --version | sed 's/bedtools //'"),
          topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id ?: tsv.simpleName}.regions"
    """
    awk -F'\\t' -v pad=${pad} -v tiers="${tiers}" '
        NR==1 {
            for (i=1; i<=NF; i++) col[\$i]=i
            n = split(tiers, tl, ",")
            for (i=1; i<=n; i++) keep[tl[i]]=1
            if (!("CHROM" in col) || !("POS" in col)) {
                print "ERROR: TSV lacks CHROM/POS columns" > "/dev/stderr"; exit 1
            }
            next
        }
        {
            t = ("driver_tier" in col) ? \$(col["driver_tier"]) : "1"
            if (!(t in keep)) next
            c = \$(col["CHROM"]); p = \$(col["POS"])
            # A breakend is a point: never span POS..END for it. See the BND note above.
            sv = ("SVTYPE" in col) ? \$(col["SVTYPE"]) : ""
            e = (sv == "BND") ? p \\
                : (("END" in col && \$(col["END"]) != "" && \$(col["END"]) != "NA") ? \$(col["END"]) : p)
            if (e < p) { tmp=p; p=e; e=tmp }        # defensive: never emit a negative interval
            s = p - 1 - pad; if (s < 0) s = 0
            print c "\\t" s "\\t" (e + pad)
        }
    ' ${tsv} | sort -k1,1 -k2,2n | bedtools merge -i - > ${prefix}.bed

    echo "INFO: \$(wc -l < ${prefix}.bed) merged intervals from \$(( \$(wc -l < ${tsv}) - 1 )) TSV rows (pad=${pad}bp, tiers=${tiers})" >&2
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id ?: tsv.simpleName}.regions"
    """
    touch ${prefix}.bed
    """
}
