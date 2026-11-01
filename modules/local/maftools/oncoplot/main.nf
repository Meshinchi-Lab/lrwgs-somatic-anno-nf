// modules/local/maftools/oncoplot/main.nf
// Render the cohort oncoprint with maftools from a MAF, a copy-number table, a gene list
// and a clinical table written by COHORT_TO_MAF. maftools has no CLI; the R program is
// written into the work directory by the script and run with Rscript (no file in bin/).
// Outputs the figure (PNG + SVG), the gene x sample matrix maftools used (onco_matrix.txt)
// and its summary tables, so the figure is reproducible from inspectable files. A MAF with
// no oncoprint-eligible rows, or a gene list with fewer than two genes, writes a NO_DATA
// marker naming the reason instead of failing (spec 8).
process MAFTOOLS_ONCOPLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bioconductor-maftools=2.26.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-maftools:2.26.0--r45h01b2380_0' :
        'quay.io/biocontainers/bioconductor-maftools:2.26.0--r45h01b2380_0' }"

    input:
    tuple val(meta), path(maf), path(cn_table), path(genes), path(clinical)

    output:
    tuple val(meta), path("${prefix}.png"),       optional: true, emit: png
    tuple val(meta), path("${prefix}.svg"),       optional: true, emit: svg
    tuple val(meta), path("onco_matrix.txt"),     optional: true, emit: matrix
    // write.mafSummary writes <prefix>_summary.txt, _geneSummary.txt, _sampleSummary.txt,
    // _clinicalData.txt and _maftools.maf: one glob, one list.
    tuple val(meta), path("${prefix}_*"),         optional: true, emit: summaries
    tuple val(meta), path("${prefix}.NO_DATA"),   optional: true, emit: no_data
    tuple val("${task.process}"), val('maftools'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"maftools\")))'"),
          topic: versions, emit: versions_maftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.oncoprint"
    """
    # The oncoprint classes: the nine classes maftools treats as non-synonymous by default
    # plus the SV classes COHORT_TO_MAF writes. Defined once, handed to awk and R as VC.
    # COHORT_TO_MAF counts the same nine classes when it builds the gene universe: keep the
    # two lists identical.
    VC='Frame_Shift_Del,Frame_Shift_Ins,Splice_Site,Translation_Start_Site,Nonsense_Mutation,Nonstop_Mutation,In_Frame_Del,In_Frame_Ins,Missense_Mutation,SV_BND,SV_DEL,SV_DUP,SV_INV,SV_INS'

    # Eligible rows = MAF data rows whose Variant_Classification is one of those classes.
    # read.maf stops with "No non-synonymous mutations found" when there are none, so a MAF
    # holding only Silent/Intron/UTR rows, or no rows at all, is NO_DATA rather than an R
    # error. A copy-number-only cohort is NO_DATA for the same reason.
    neligible=\$(awk -F'\\t' -v vc="\${VC}" '
        BEGIN { n = split(vc, a, ","); for (i = 1; i <= n; i++) ok[a[i]] = 1 }
        /^#/ { next }
        !hdr { for (i = 1; i <= NF; i++) if (\$i == "Variant_Classification") c = i; hdr = 1; next }
        c && (\$c in ok) { k++ }
        END {
            if (!c) { print "ERROR: MAF has no Variant_Classification column" > "/dev/stderr"; exit 1 }
            print k + 0
        }' ${maf})
    # The gene list is the COHORT_TO_MAF row rule applied as-is (spec 4/5): fewer than two
    # genes is NO_DATA (maftools oncoplot needs two rows), never a frequency-ranked fallback.
    # Non-blank lines are counted the way R trims them below. Read through redirection so a
    # missing file, or [] (empty name), fails here instead of busybox awk reading stdin.
    ngenes=\$(awk 'NF { n++ } END { print n + 0 }' < "${genes}")
    # cn_table and clinical may be [] (nothing staged, empty name) or header-only; both count
    # as absent and are not passed to read.maf. awk NR counts an unterminated last line,
    # which wc -l would miss.
    ncn=0;   if [ -s "${cn_table}" ]; then ncn=\$(awk 'END { print NR }' "${cn_table}"); fi
    nclin=0; if [ -s "${clinical}" ]; then nclin=\$(awk 'END { print NR }' "${clinical}"); fi
    echo "INFO: eligible MAF rows=\${neligible} genes=\${ngenes} cn_lines=\${ncn} clinical_lines=\${nclin}" >&2

    if [ "\${neligible}" -eq 0 ]; then
        echo "no oncoprint-eligible MAF rows (classes: \${VC})" > ${prefix}.NO_DATA
    elif [ "\${ngenes}" -eq 0 ]; then
        echo "no oncoprint-eligible genes (the gene list is empty)" > ${prefix}.NO_DATA
    elif [ "\${ngenes}" -lt 2 ]; then
        echo "only one oncoprint gene; maftools oncoplot needs at least two" > ${prefix}.NO_DATA
    else
        # The R program is static text (quoted heredoc: no shell expansion, no quoting rules)
        # and takes every file name and count from the environment. It stays in the work
        # directory, so a failed task can be re-run by hand with the same inputs. The EOF
        # terminator sits at the script's base indentation: Nextflow strips that common
        # indent, which puts EOF at column 0 where bash requires it.
        cat > oncoplot.R <<'EOF'
        suppressPackageStartupMessages(library(maftools))
        e  <- Sys.getenv(c("VC", "MAF", "CN", "CLIN", "GENES", "PREFIX", "NCN", "NCLIN"))
        vc <- strsplit(e[["VC"]], ",")[[1]]
        # Optional tables: absent or header-only means not passed. The clinical table is read
        # as character so numerically coded sex/etp (0/1, 1/2) stay categorical; maftools
        # would otherwise treat them as continuous and demand a sequential palette.
        cn   <- if (as.integer(e[["NCN"]])   > 1) e[["CN"]] else NULL
        clin <- if (as.integer(e[["NCLIN"]]) > 1)
                  data.table::fread(e[["CLIN"]], sep = "\\t", colClasses = "character") else NULL
        # removeDuplicatedVariants = FALSE is REQUIRED: maftools de-duplicates on
        # (Chromosome, Start, Sample, Ref, Alt) without Hugo_Symbol, which would silently drop
        # every gene after the first on a multi-gene SV row (see the COHORT_TO_MAF header).
        m <- read.maf(e[["MAF"]], cnTable = cn, clinicalData = clin, vc_nonSyn = vc,
                      removeDuplicatedVariants = FALSE, verbose = FALSE)
        genes <- trimws(readLines(e[["GENES"]])); genes <- genes[nzchar(genes)]
        # Gene names sit in the left margin, which maftools fixes at 5 lines: long symbols
        # (LINC00649, TMPRSS15) were clipped. Size the margin from the longest name instead.
        gmar  <- max(5, ceiling(0.8 * max(nchar(genes))) + 1)
        # The maftools 2.26.0 defaults for its own classes, written out because
        # get_vcColors() is not exported, so this figure matches every other maftools
        # figure. The SV classes take hues the defaults leave free.
        cols <- c(Frame_Shift_Del = "#1F78B4", Frame_Shift_Ins = "#6A3D9A", Splice_Site = "#FF7F00",
                  Translation_Start_Site = "#F46D43", Nonsense_Mutation = "#E31A1C", Nonstop_Mutation = "#A6CEE3",
                  In_Frame_Del = "#FFFF99", In_Frame_Ins = "#D53E4F", Missense_Mutation = "#33A02C",
                  Multi_Hit = "#000000", Amp = "#EE82EE", Del = "#4169E1",
                  # The bar charts count a gene with both a copy-number change and an SNV/SV
                  # alteration in one sample as Complex_Event (maftools default colour).
                  Complex_Event = "#7B7060",
                  SV_BND = "#E7298A", SV_DEL = "#1B9E77", SV_DUP = "#BCBD22", SV_INV = "#17BECF", SV_INS = "#8C564B")
        # Clinical tracks: one palette per feature (warm for sex, cool for etp, no colour in
        # both), assigned to the sorted (C locale) non-empty observed values; empty strings
        # and NA are left to maftools, which draws them as missing, and a feature with no
        # non-empty value is not drawn. More levels than colours is a data error (free text
        # in sex or etp) and is reported rather than recycled into a misleading legend.
        # maftools 2.26.0 cannot draw clinical tracks for a one-sample oncoprint (its
        # annotation matrix collapses to a vector), so a cohort of one gets none.
        pal   <- list(sex = c("#FC8D62", "#FFD92F", "#E5C494", "#B15928", "#F781BF", "#FDB462"),
                      etp = c("#8DA0CB", "#A6D854", "#66C2A5", "#80B1D3", "#B3B3B3", "#984EA3"))
        feats <- character(0); acol <- list()
        if (!is.null(clin) && nrow(getClinicalData(m)) >= 2) {
          cd <- getClinicalData(m)
          for (f in intersect(names(pal), names(cd))) {
            v <- as.character(cd[[f]]); v <- sort(unique(v[!is.na(v) & nzchar(v)]), method = "radix")
            if (length(v) == 0) next
            if (length(v) > length(pal[[f]]))
              stop("clinical feature ", f, " has ", length(v), " levels (", paste(v, collapse = ","),
                   "); at most ", length(pal[[f]]), " are supported")
            feats <- c(feats, f); acol[[f]] <- setNames(pal[[f]][seq_along(v)], v)
          }
        }
        if (length(feats) == 0) { feats <- NULL; acol <- NULL }
        draw <- function() oncoplot(m, genes = genes, top = length(genes), colors = cols,
                                    clinicalFeatures = feats, annotationColor = acol,
                                    showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE,
                                    drawRowBar = TRUE, drawColBar = TRUE, writeMatrix = TRUE,
                                    gene_mar = gmar)
        png(paste0(e[["PREFIX"]], ".png"), width = 2400, height = 1600, res = 200); draw(); invisible(dev.off())
        svg(paste0(e[["PREFIX"]], ".svg"), width = 12, height = 8); draw(); invisible(dev.off())
        write.mafSummary(m, basename = e[["PREFIX"]])
    EOF
        VC="\${VC}" MAF="${maf}" CN="${cn_table}" CLIN="${clinical}" GENES="${genes}" \\
            PREFIX="${prefix}" NCN="\${ncn}" NCLIN="\${nclin}" Rscript oncoplot.R
        # Every output is optional (the NO_DATA path), so a missing figure would not fail the
        # task by itself: check explicitly. R errors already fail it through set -e.
        for f in ${prefix}.png ${prefix}.svg onco_matrix.txt; do
            [ -s "\$f" ] || { echo "ERROR: maftools wrote no \$f" >&2; exit 1; }
        done
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.oncoprint"
    """
    printf 'stub' > ${prefix}.png
    printf '<svg xmlns="http://www.w3.org/2000/svg"></svg>' > ${prefix}.svg
    printf '\\tstub\\n' > onco_matrix.txt
    for s in summary geneSummary sampleSummary clinicalData; do : > ${prefix}_\${s}.txt; done
    : > ${prefix}_maftools.maf
    """
}
