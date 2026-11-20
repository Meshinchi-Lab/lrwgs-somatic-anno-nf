// modules/local/cohort_to_maf/main.nf
// Reshape the pass-1 cohort TSVs into maftools inputs: a MAF (SNV + SV rows), a
// copy-number table, an optional clinical table and the oncoprint gene list.
// csvtk selects columns by name; awk does the reshaping. No Python/R.
//
// Rules (spec 6.2):
//   SNV rows : SYMBOL present; class from the VEP lookup (first `&` term); sentinels
//              Frame_Shift_BY_TYPE / In_Frame_BY_LENGTH resolved from allele lengths;
//              Variant_Type SNP/DNP/TNP/ONP/INS/DEL from lengths; unmapped -> sidecar.
//              Indels follow the MAF convention: the shared anchor base is dropped, a
//              deletion starts one base after the VCF POS with Tumor_Seq_Allele2 "-", an
//              insertion has Reference_Allele "-" and End = Start + 1.
//   SV rows  : class SV_<SVTYPE> (upper-cased, TRA -> BND); Variant_Type DEL/INS/SV.
//              BND -> the left breakend gene at the record's own locus, with its REF/ALT.
//              A BND pair is two records (SAMPLE + ID minus _1/_2; the SV projection is read
//              twice to pair them). With both present each writes its own row, and the _2 row
//              is dropped when its gene equals the _1 record's (intragenic: one row). A record
//              whose mate is absent -- AnnotSV keeps one record of a same-chromosome pair --
//              also writes the right breakend gene AT THE PARTNER LOCUS (from `breakpoint`)
//              with REF "N" and the mate's ALT, mirrored to point back at the record:
//              t[p[ -> ]c:p]N, t]p] -> N]c:p], ]p]t -> N[c:p[, [p[t -> [c:p[N. Genes only
//              with gene_source == overlap. The rows do not depend on which mate comes first.
//              Other types -> each gene in Candidate_genes (panel-intersected upstream),
//              deduplicated within the record.
//   CN table : Gene, Sample_name, Amp|Del for Candidate_genes of each CNA segment; a gene
//              with both states in one sample keeps the first and is counted (cn_conflicts).
//   genes    : (panel & altered) U (altered in >= 2 samples), by altered-sample count
//              desc, panel first, alphabetical (C locale); capped at top_genes. No panel ->
//              all altered genes by frequency. "Altered" counts only the classes maftools
//              treats as non-synonymous, the SV_* classes and copy number.
//   rows with an empty SAMPLE or a non-numeric coordinate are skipped and counted
//   (skipped_no_sample, skipped_bad_coordinate); a MAF row with a non-numeric position
//   would break read.maf. skipped_sv_duplicate_mates counts a repeated SAMPLE + ID and the
//   _2 row of an intragenic BND pair.
//
// Consumer contract (MAFTOOLS_ONCOPLOT): read.maf must be called with
// removeDuplicatedVariants = FALSE -- maftools de-duplicates on
// (Chromosome, Start, Sample, Ref, Alt) WITHOUT Hugo_Symbol, which would silently drop
// every gene after the first on a multi-gene DEL/DUP/INV row. Two or more non-CN rows for
// one gene in one sample render as Multi_Hit; that is maftools' convention.
//
// Limits: busybox awk/sort; csvtk re-quotes a field containing a double quote (the report
// tables carry none); the clinical table is read from the SV table only, so a sample with
// CNA/SNV rows but no SV row gets no clinical row.
process COHORT_TO_MAF {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::csvtk=0.22.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.22.0--h9ee0642_1' :
        'quay.io/biocontainers/csvtk:0.22.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(sv_tsv)
    path cna_tsv            // or []
    path snv_tsv            // or []
    path candidate_genes    // one symbol per line, or []
    path vep_lookup         // assets/genomic_viz/vep_consequence_to_maf.tsv
    val  top_genes

    output:
    tuple val(meta), path("cohort.maf"),                 emit: maf
    tuple val(meta), path("cohort.cn.tsv"),              emit: cn_table
    tuple val(meta), path("cohort.clinical.tsv"),        emit: clinical
    tuple val(meta), path("oncoprint_genes.txt"),        emit: genes
    tuple val(meta), path("unmapped_consequences.tsv"),  emit: unmapped
    tuple val(meta), path("maf.counts.tsv"),             emit: counts
    tuple val("${task.process}"), val('csvtk'), eval("csvtk version | sed 's/csvtk v//'"),
          topic: versions, emit: versions_csvtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def cna_in   = cna_tsv         ? cna_tsv.toString()         : ''
    def snv_in   = snv_tsv         ? snv_tsv.toString()         : ''
    def panel_in = candidate_genes ? candidate_genes.toString() : ''
    def top      = top_genes.toString()
    if ( !(top ==~ /[1-9][0-9]*/) ) {
        error("COHORT_TO_MAF: top_genes must be a positive integer, got '${top}'")
    }
    """
    set -euo pipefail

    # One csvtk call per file, matched with a here-string: no pipe, so grep -q's early
    # exit cannot SIGPIPE csvtk and turn pipefail into a false "missing column".
    require_cols() {
        f="\$1"; shift
        hdr="\$(csvtk headers -t "\$f")" || { echo "ERROR: cannot read a header from \$f" >&2; exit 1; }
        missing=""
        for c in "\$@"; do
            grep -qxF -- "\$c" <<< "\$hdr" || missing="\$missing \$c"
        done
        if [ -n "\$missing" ]; then
            echo "ERROR: \$f lacks required column(s):\$missing" >&2
            exit 1
        fi
    }
    require_cols ${sv_tsv} SAMPLE ID CHROM POS END REF ALT SVTYPE breakpoint driver_tier Candidate_genes gene_name_left gene_name_right gene_source_left gene_source_right
    [ -z "${cna_in}" ] || require_cols ${cna_in} SAMPLE CHROM POS END WAKHAN_TYPE driver_tier Candidate_genes
    [ -z "${snv_in}" ] || require_cols ${snv_in} SAMPLE CHROM POS REF ALT SYMBOL Consequence driver_tier

    # ---- projections --------------------------------------------------------------
    csvtk cut -t -f SAMPLE,ID,CHROM,POS,END,REF,ALT,SVTYPE,driver_tier,Candidate_genes,gene_name_left,gene_name_right,gene_source_left,gene_source_right,breakpoint ${sv_tsv} > sv.proj.tsv
    if [ -n "${cna_in}" ]; then
        csvtk cut -t -f SAMPLE,CHROM,POS,END,WAKHAN_TYPE,driver_tier,Candidate_genes ${cna_in} > cna.proj.tsv
    else
        printf 'SAMPLE\\tCHROM\\tPOS\\tEND\\tWAKHAN_TYPE\\tdriver_tier\\tCandidate_genes\\n' > cna.proj.tsv
    fi
    if [ -n "${snv_in}" ]; then
        csvtk cut -t -f SAMPLE,CHROM,POS,REF,ALT,SYMBOL,Consequence,driver_tier ${snv_in} > snv.proj.tsv
    else
        printf 'SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tSYMBOL\\tConsequence\\tdriver_tier\\n' > snv.proj.tsv
    fi
    # Clinical columns are optional; the SV table carries them when sample metadata was supplied.
    svhdr="\$(csvtk headers -t ${sv_tsv})"
    if grep -qxF sex <<< "\$svhdr" && grep -qxF etp <<< "\$svhdr"; then
        csvtk cut -t -f SAMPLE,sex,etp ${sv_tsv} > clin.proj.tsv
    else
        printf 'SAMPLE\\tsex\\tetp\\n' > clin.proj.tsv
    fi
    if [ -n "${panel_in}" ]; then grep -v '^#' ${panel_in} | tr -d '\\r' | sed '/^\\s*\$/d' > panel.txt; else : > panel.txt; fi

    # ---- MAF + CN + clinical + gene/sample pairs ----------------------------------
    awk -F'\\t' -v OFS='\\t' \\
        -v lookf=${vep_lookup} -v svf=sv.proj.tsv -v cnaf=cna.proj.tsv -v snvf=snv.proj.tsv -v clinf=clin.proj.tsv '
        function chrom(c) { if (c !~ /^chr/) c = "chr" c; return c }
        function isnum(x){ return (x ~ /^[0-9]+\$/) }
        function vtype(r, a) {
            lr = length(r); la = length(a)
            if (lr == la) { if (lr == 1) return "SNP"; if (lr == 2) return "DNP"; if (lr == 3) return "TNP"; return "ONP" }
            return (la > lr) ? "INS" : "DEL"
        }
        # "altered" for the gene universe: the classes maftools treats as non-synonymous
        # (its vc_nonSyn default), every SV_* class, and copy number (added by the CNA pass).
        function nonsyn(vc) { return (vc in NS) || (vc ~ /^SV_/) }
        function maf_row(gene, chr, s, e, ref, alt, vc, vt, sample, tier, src) {
            print gene, chrom(chr), s, e, ref, alt, vc, vt, sample, tier, src > "cohort.maf"
            if (nonsyn(vc)) print gene, sample > "pairs.tsv"
            nmaf[src]++; nclass[vc]++
        }
        function trim(s) { gsub(/^ +| +\$/, "", s); return s }
        function lgene() { return (\$13 == "overlap" && \$11 != "" && \$11 != "NA") ? \$11 : "" }
        function rgene() { return (\$14 == "overlap" && \$12 != "" && \$12 != "NA") ? \$12 : "" }
        # mate number from the _1 / _2 suffix of a BND ID ("" when it has none)
        function mate(id) { return (id ~ /_1\$/) ? 1 : (id ~ /_2\$/) ? 2 : "" }
        # ALT of the partner breakend, pointing back at the locus cp of this record (VCF 4.2, 5.4):
        # t[p[ -> ]cp]N, t]p] -> N]cp], ]p]t -> N[cp[, [p[t -> [cp[N; a symbolic ALT is kept.
        function mate_alt(alt, cp,    f, l) {
            f = substr(alt, 1, 1); l = substr(alt, length(alt), 1)
            if (l == "[") return "]" cp "]N"
            if (l == "]") return "N]" cp "]"
            if (f == "]") return "N[" cp "["
            if (f == "[") return "[" cp "[N"
            return alt
        }
        BEGIN {
            split("Frame_Shift_Del Frame_Shift_Ins Splice_Site Translation_Start_Site Nonsense_Mutation Nonstop_Mutation In_Frame_Del In_Frame_Ins Missense_Mutation", nsl, " ")
            for (q in nsl) NS[nsl[q]] = 1
            print "#version 2.4" > "cohort.maf"
            print "Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","driver_tier","Alteration_Source" > "cohort.maf"
            print "Gene","Sample_name","CN" > "cohort.cn.tsv"
            print "Tumor_Sample_Barcode","sex","etp" > "cohort.clinical.tsv"
        }
        FILENAME == lookf { if (FNR > 1) cls[\$1] = \$2; next }

        # SNV: 1 SAMPLE 2 CHROM 3 POS 4 REF 5 ALT 6 SYMBOL 7 Consequence 8 tier
        FILENAME == snvf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            if (\$6 == "" || \$6 == "NA" || \$6 == ".") { skip_sym++; next }
            if (!isnum(\$3)) { badcoord++; next }
            split(\$7, terms, "&"); term = terms[1]
            if (!(term in cls)) { unm[(term == "") ? "(empty)" : term]++; next }
            r = \$4; a = \$5; vc = cls[term]; vt = vtype(r, a); d = length(a) - length(r)
            if (vc == "Frame_Shift_BY_TYPE") vc = (vt == "INS") ? "Frame_Shift_Ins" : "Frame_Shift_Del"
            if (vc == "In_Frame_BY_LENGTH") {
                if (d == 0) vc = "Missense_Mutation"; else if (d % 3 == 0) vc = (d > 0) ? "In_Frame_Ins" : "In_Frame_Del"; else vc = "Missense_Mutation"
            }
            s = \$3 + 0
            if (vt == "INS" || vt == "DEL") {
                # MAF convention: drop the shared anchor base (if any) and shift the start past it
                tr = 0
                if (substr(r, 1, 1) == substr(a, 1, 1)) { r = substr(r, 2); a = substr(a, 2); tr = 1 }
                if (r == "") r = "-"
                if (a == "") a = "-"
                if (r == "-") e = s + 1                       # pure insertion: between the flanking bases
                else { if (tr) s = s + 1; e = s + length(r) - 1 }
            } else { e = s + length(r) - 1 }
            maf_row(\$6, \$2, s, e, r, a, vc, vt, \$1, \$8, "SNV")
            next
        }

        # SV: 1 SAMPLE 2 ID 3 CHROM 4 POS 5 END 6 REF 7 ALT 8 SVTYPE 9 tier 10 Candidate_genes
        #     11 gene_left 12 gene_right 13 src_left 14 src_right 15 breakpoint
        # Read twice: the first read records which BND mates are valid, and the _1 left gene.
        FILENAME == svf && FNR == 1 { svpass++; next }
        FILENAME == svf && svpass == 1 {
            st = toupper(\$8)
            if (\$1 == "" || (st != "BND" && st != "TRA") || !isnum(\$4)) next
            m = mate(\$2); if (m == "") next
            stem = \$2; sub(/_[12]\$/, "", stem)
            valid[\$1, stem, m] = 1
            if (m == 1) lg1[\$1, stem] = lgene()
            next
        }
        FILENAME == svf {
            if (\$1 == "") { nosample++; next }
            st = toupper(\$8); if (st == "TRA") st = "BND"
            idk = \$1 SUBSEP \$2
            if (idk in seen) { dup++; next }
            if (st != "BND" && st != "DEL" && st != "DUP" && st != "INV" && st != "INS") { badsv++; next }
            if (!isnum(\$4) || (st != "BND" && !isnum(\$5))) { badcoord++; next }
            seen[idk] = 1                              # only a VALID record blocks a repeat
            vc = "SV_" st
            vt = (st == "DEL") ? "DEL" : (st == "INS") ? "INS" : "SV"
            ref = (\$6 == "" || \$6 == ".") ? "N" : \$6
            if (st == "BND") {
                lg = lgene(); rg = rgene()
                m = mate(\$2); stem = \$2; sub(/_[12]\$/, "", stem)
                paired = (m != "" && ((\$1, stem, 3 - m) in valid))
                # intragenic pair: the _1 row already names the gene
                if (paired && m == 2 && lg != "" && lg == lg1[\$1, stem]) { dup++; next }
                if (lg != "") maf_row(lg, \$3, \$4, \$4, ref, \$7, vc, vt, \$1, \$9, "SV")
                if (!paired && rg != "" && rg != lg) {
                    # the partner locus: breakpoint is chrA:posA|chrB:posB
                    pc = ""; pp = ""
                    if (split(\$15, bp, "|") == 2) { split(bp[2], b, ":"); pc = b[1]; pp = b[2] }
                    if (pc == "" || !isnum(pp)) badcoord++
                    else maf_row(rg, pc, pp, pp, "N", mate_alt(\$7, \$3 ":" \$4), vc, vt, \$1, \$9, "SV")
                }
            } else if (\$10 != "" && \$10 != "NA") {
                delete gseen
                n = split(\$10, g, ";")
                for (i = 1; i <= n; i++) {
                    gg = trim(g[i]); if (gg == "" || (gg in gseen)) continue
                    gseen[gg] = 1
                    maf_row(gg, \$3, \$4, \$5, ref, \$7, vc, vt, \$1, \$9, "SV")
                }
            }
            next
        }

        # CNA: 1 SAMPLE 2 CHROM 3 POS 4 END 5 WAKHAN_TYPE 6 tier 7 Candidate_genes
        FILENAME == cnaf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            ty = toupper(\$5)
            if (ty == "GAIN" || ty == "AMP") cn = "Amp"; else if (ty == "LOSS" || ty == "DEL") cn = "Del"; else { badcna++; next }
            if (\$7 == "" || \$7 == "NA") next
            n = split(\$7, g, ";")
            for (i = 1; i <= n; i++) {
                gg = trim(g[i]); if (gg == "") continue
                k = gg SUBSEP \$1
                if (k in cnstate) { if (cnstate[k] != cn) cnconf++; continue }
                cnstate[k] = cn
                print gg, \$1, cn > "cohort.cn.tsv"; print gg, \$1 > "pairs.tsv"; ncn++
            }
            next
        }

        # clinical: 1 SAMPLE 2 sex 3 etp -- first informative row per sample; a literal NA
        # (readr output for a sample missing from the metadata) counts as empty
        FILENAME == clinf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            if (\$1 in clseen) next
            sx = (\$2 == "NA") ? "" : \$2; et = (\$3 == "NA") ? "" : \$3
            if (sx == "" && et == "") next
            clseen[\$1] = 1
            print \$1, sx, et > "cohort.clinical.tsv"; nclin++
            next
        }

        END {
            for (t in unm) print t, unm[t] > "unmapped.tmp"
            print "metric", "value" > "maf.counts.tsv"
            print "maf_rows_snv", nmaf["SNV"] + 0 > "maf.counts.tsv"
            print "maf_rows_sv", nmaf["SV"] + 0 > "maf.counts.tsv"
            print "cn_rows", ncn + 0 > "maf.counts.tsv"
            print "cn_conflicts", cnconf + 0 > "maf.counts.tsv"
            print "clinical_rows", nclin + 0 > "maf.counts.tsv"
            print "skipped_snv_no_symbol", skip_sym + 0 > "maf.counts.tsv"
            print "skipped_sv_duplicate_mates", dup + 0 > "maf.counts.tsv"
            print "skipped_sv_unknown_type", badsv + 0 > "maf.counts.tsv"
            print "skipped_cna_unknown_type", badcna + 0 > "maf.counts.tsv"
            print "skipped_bad_coordinate", badcoord + 0 > "maf.counts.tsv"
            print "skipped_no_sample", nosample + 0 > "maf.counts.tsv"
            nu = 0; for (t in unm) nu++
            print "unmapped_consequence_terms", nu > "maf.counts.tsv"
            for (v in nclass) print "class_" v, nclass[v] > "class.tmp"
            printf("INFO: MAF rows SNV=%d SV=%d; CN rows=%d; clinical=%d; unmapped terms=%d\\n", nmaf["SNV"] + 0, nmaf["SV"] + 0, ncn + 0, nclin + 0, nu) > "/dev/stderr"
            if (nosample + badcoord > 0) printf("WARN: skipped %d rows with an empty SAMPLE and %d with a non-numeric coordinate\\n", nosample + 0, badcoord + 0) > "/dev/stderr"
        }' ${vep_lookup} snv.proj.tsv sv.proj.tsv sv.proj.tsv cna.proj.tsv clin.proj.tsv
    [ -f pairs.tsv ] || : > pairs.tsv
    [ -f unmapped.tmp ] || : > unmapped.tmp
    [ -f class.tmp ] || : > class.tmp
    # awk hash order is not stable: sort the per-class and unmapped rows (C locale).
    { printf 'vep_term\\tn_rows\\n'; LC_ALL=C sort -k1,1 unmapped.tmp; } > unmapped_consequences.tsv
    LC_ALL=C sort -k1,1 class.tmp >> maf.counts.tsv

    # ---- gene universe ------------------------------------------------------------
    # distinct (gene, sample) pairs -> per-gene sample count; panel membership as a sort key.
    # The cap is applied by awk, not head: head exits early and SIGPIPEs the sort under pipefail.
    LC_ALL=C sort -u pairs.tsv | awk -F'\\t' -v OFS='\\t' -v panelf=panel.txt '
        FILENAME == panelf { panel[\$1] = 1; np++; next }
        { n[\$1]++ }
        END {
            for (g in n) {
                inpanel = (g in panel) ? 1 : 0
                keep = (np == 0) ? 1 : (inpanel || n[g] >= 2)
                if (keep) print g, n[g], inpanel
            }
        }' panel.txt - \\
        | LC_ALL=C sort -t "\$(printf '\\t')" -k2,2nr -k3,3nr -k1,1 > universe.tsv
    awk -F'\\t' -v n=${top} 'NR <= n { print \$1 }' universe.tsv > oncoprint_genes.txt
    # The report caption states how many genes qualified and the cap that was applied.
    # Assigned first so a failing awk stops the task under set -e (inside printf's arguments it would not).
    n_eligible=\$(awk 'END { print NR }' universe.tsv)
    printf 'oncoprint_genes_eligible\\t%s\\noncoprint_genes_cap\\t%s\\n' "\${n_eligible}" "${top}" >> maf.counts.tsv
    rm -f sv.proj.tsv cna.proj.tsv snv.proj.tsv clin.proj.tsv panel.txt pairs.tsv unmapped.tmp class.tmp universe.tsv
    """

    stub:
    """
    printf '#version 2.4\\nHugo_Symbol\\tChromosome\\tStart_Position\\tEnd_Position\\tReference_Allele\\tTumor_Seq_Allele2\\tVariant_Classification\\tVariant_Type\\tTumor_Sample_Barcode\\tdriver_tier\\tAlteration_Source\\n' > cohort.maf
    printf 'Gene\\tSample_name\\tCN\\n' > cohort.cn.tsv
    printf 'Tumor_Sample_Barcode\\tsex\\tetp\\n' > cohort.clinical.tsv
    : > oncoprint_genes.txt
    printf 'vep_term\\tn_rows\\n' > unmapped_consequences.tsv
    printf 'metric\\tvalue\\n' > maf.counts.tsv
    """
}
