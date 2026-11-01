// modules/local/cohort_to_circos/main.nf
// Reshape the report's pass-1 cohort TSVs into Circos text tracks: one set for the
// cohort and one per sample. csvtk selects columns BY NAME (the real TSVs have ~100
// columns; the fixtures have 18), awk does the relational work. No Python/R.
//
// Encodings (spec 6.1):
//   links   : one line per SV; BND mates collapsed on SAMPLE + ID minus _1/_2; INS skipped;
//             colour = tier<k>; thickness from n_any (1->2p, 2->4p, 3-4->6p, >=5->8p;
//             SV absent from the map, or no map -> 3p); z = 4 - tier so Tier 1 draws on top.
//             for a BND the larger n_any of its two mates is used (Jasmine clusters each
//             mate separately), so thickness does not depend on which mate is listed first.
//   cna     : tiles coloured gain|loss from WAKHAN_TYPE (GAIN/AMP, LOSS/DEL).
//   labels  : Tier-1 gene names -- SV breakend genes with gene_source == overlap
//             (fallback: first Gene_name_display token), CNA Candidate_genes at the
//             segment midpoint, SNV SYMBOL. Deduplicated per set on (gene, chr, 1 Mb bin).
//   contigs : chr prefix stripped, hs prefix added, only 1-22/X/Y kept.
//   tiers   : driver_tier outside 1..3 (e.g. NA) is drawn as tier 3.
//   rows with a non-numeric coordinate count as bad_breakpoint; rows with an empty SAMPLE
//   are skipped and reported on stderr.
//
// Limits: busybox awk keeps every output file open (3 per sample + 4); at ulimit -n 1024
// that is ~340 samples, failing loudly with "Too many open files". csvtk re-quotes a field
// containing a double quote; the report tables carry none.
process COHORT_TO_CIRCOS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::csvtk=0.22.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.22.0--h9ee0642_1' :
        'quay.io/biocontainers/csvtk:0.22.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(sv_tsv)
    path cna_tsv            // cohort_cna_pass2.tsv or []
    path snv_tsv            // cohort_snv.tsv or []
    path recurrence_map     // cohort_sv_recurrent.map.tsv or []

    output:
    tuple val(meta), path("*.links.txt"),  emit: links
    tuple val(meta), path("*.cna.txt"),    emit: cna
    tuple val(meta), path("*.labels.txt"), emit: labels
    tuple val(meta), path("counts.tsv"),   emit: counts
    tuple val("${task.process}"), val('csvtk'), eval("csvtk version | sed 's/csvtk v//'"),
          topic: versions, emit: versions_csvtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def cna_in = cna_tsv        ? cna_tsv.toString()        : ''
    def snv_in = snv_tsv        ? snv_tsv.toString()        : ''
    def map_in = recurrence_map ? recurrence_map.toString() : ''
    """
    set -euo pipefail

    # ---- 1. Required columns: fail fast, naming every missing column -----------------
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
    require_cols ${sv_tsv} SAMPLE ID CHROM POS END SVTYPE breakpoint driver_tier Gene_name_display gene_name_left gene_name_right gene_source_left gene_source_right
    [ -z "${cna_in}" ] || require_cols ${cna_in} SAMPLE CHROM POS END WAKHAN_TYPE driver_tier Candidate_genes
    [ -z "${snv_in}" ] || require_cols ${snv_in} SAMPLE CHROM POS SYMBOL driver_tier

    # ---- 2. Project to fixed column orders so awk can use positional fields -----------
    csvtk cut -t -f SAMPLE,ID,CHROM,POS,END,SVTYPE,breakpoint,driver_tier,Gene_name_display,gene_name_left,gene_name_right,gene_source_left,gene_source_right ${sv_tsv} > sv.proj.tsv
    if [ -n "${cna_in}" ]; then
        csvtk cut -t -f SAMPLE,CHROM,POS,END,WAKHAN_TYPE,driver_tier,Candidate_genes ${cna_in} > cna.proj.tsv
    else
        printf 'SAMPLE\\tCHROM\\tPOS\\tEND\\tWAKHAN_TYPE\\tdriver_tier\\tCandidate_genes\\n' > cna.proj.tsv
    fi
    if [ -n "${snv_in}" ]; then
        csvtk cut -t -f SAMPLE,CHROM,POS,SYMBOL,driver_tier ${snv_in} > snv.proj.tsv
    else
        printf 'SAMPLE\\tCHROM\\tPOS\\tSYMBOL\\tdriver_tier\\n' > snv.proj.tsv
    fi
    if [ -n "${map_in}" ]; then
        csvtk cut -t -f source_id,n_any ${map_in} > rec.proj.tsv
    else
        printf 'source_id\\tn_any\\n' > rec.proj.tsv
    fi

    # ---- 3. Sample set: every SAMPLE seen in any table gets a (possibly empty) file set ---
    # tr's class must equal awk's san() below, or a sample's files would be written under a
    # different name than the pre-created empty ones.
    awk -F'\\t' 'FNR>1 && \$1 != "" { print \$1 }' sv.proj.tsv cna.proj.tsv snv.proj.tsv | sort -u > samples.raw.txt
    tr -c 'A-Za-z0-9._\\n-' '_' < samples.raw.txt | sort -u > samples.txt
    if grep -qx cohort samples.txt; then
        echo "ERROR: a SAMPLE named 'cohort' collides with the cohort track set; rename the sample" >&2
        exit 1
    fi
    nraw=\$(wc -l < samples.raw.txt); nsan=\$(wc -l < samples.txt)
    [ "\$nraw" -eq "\$nsan" ] || echo "WARN: \$nraw distinct SAMPLE values collapsed to \$nsan sanitised set names" >&2
    for s in cohort \$(cat samples.txt); do
        : > "\${s}.links.txt"; : > "\${s}.cna.txt"; : > "\${s}.labels.txt"
    done

    # ---- 4. One awk pass over the four projections ----------------------------------
    awk -F'\\t' -v OFS='\\t' \\
        -v recf=rec.proj.tsv -v svf=sv.proj.tsv -v cnaf=cna.proj.tsv -v snvf=snv.proj.tsv '
        function san(s)  { gsub(/[^A-Za-z0-9._-]/, "_", s); return s }
        function norm(c) { sub(/^chr/, "", c); return c }
        function ok(c)   { return (c ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)\$/) }
        function isnum(x){ return (x ~ /^[0-9]+\$/) }
        function thick(n){ if (n == "") return 3; n = n + 0; if (n <= 1) return 2; if (n == 2) return 4; if (n <= 4) return 6; return 8 }
        function tierof(t){ t = t + 0; if (t < 1 || t > 3) return 3; return t }
        function label(set, chr, pos, gene) {
            if (gene == "" || gene == "NA" || gene == ".") return
            k = set SUBSEP gene SUBSEP chr SUBSEP int(pos / 1000000)
            if (k in lseen) return
            lseen[k] = 1
            print "hs" chr, pos, pos, gene > (set ".labels.txt")
            nlab[set]++
        }
        function emit_labels(chr, pos, gene) { label("cohort", chr, pos, gene); label(ss, chr, pos, gene) }
        BEGIN { OFS = " " }

        # recurrence map: source_id -> n_any, plus the per-mate-pair maximum
        FILENAME == recf {
            if (FNR > 1) {
                nany[\$1] = \$2
                stem = \$1; sub(/_[12]\$/, "", stem)
                if (!(stem in nmax) || \$2 + 0 > nmax[stem]) nmax[stem] = \$2 + 0
            }
            next
        }

        # SV projection: 1 SAMPLE 2 ID 3 CHROM 4 POS 5 END 6 SVTYPE 7 breakpoint 8 tier
        #                9 Gene_name_display 10 gene_left 11 gene_right 12 src_left 13 src_right
        FILENAME == svf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            ss = san(\$1); sets[ss] = 1
            st = toupper(\$6); if (st == "TRA") st = "BND"      # same normalisation as COHORT_TO_MAF
            if (st == "INS") { skip_ins[ss]++; next }
            key = \$2; if (st == "BND") sub(/_[12]\$/, "", key); key = ss SUBSEP key
            if (key in seen) { dup[ss]++; next }
            if (st == "BND") {
                n = split(\$7, bp, "|")
                if (n != 2) { badbp[ss]++; next }
                split(bp[1], a, ":"); split(bp[2], b, ":")
                c1 = norm(a[1]); p1 = a[2]; c2 = norm(b[1]); p2 = b[2]
            } else {
                c1 = norm(\$3); p1 = \$4; c2 = c1; p2 = \$5
            }
            if (!isnum(p1) || !isnum(p2)) { badbp[ss]++; next }
            p1 = p1 + 0; p2 = p2 + 0
            if (st != "BND" && p2 < p1) { t = p1; p1 = p2; p2 = t }
            if (!ok(c1) || !ok(c2)) { offkar[ss]++; next }
            seen[key] = 1                              # only a DRAWN record blocks its mate
            tier = tierof(\$8)
            rec = \$1 "_" \$2; stem = rec; if (st == "BND") sub(/_[12]\$/, "", stem)
            w = thick((st == "BND" && (stem in nmax)) ? nmax[stem] : ((rec in nany) ? nany[rec] : ""))
            line = sprintf("hs%s %d %d hs%s %d %d color=tier%d,thickness=%dp,z=%d", c1, p1, p1, c2, p2, p2, tier, w, 4 - tier)
            print line > "cohort.links.txt"; print line > (ss ".links.txt")
            nlinks["cohort"]++; nlinks[ss]++
            if (tier == 1) {
                got = 0
                if (\$12 == "overlap" && \$10 != "" && \$10 != "NA") { emit_labels(c1, p1, \$10); got = 1 }
                if (\$13 == "overlap" && \$11 != "" && \$11 != "NA") { emit_labels(c2, p2, \$11); got = 1 }
                if (!got) { split(\$9, g, ";"); gsub(/^ +| +\$/, "", g[1]); emit_labels(c1, p1, g[1]) }
            }
            next
        }

        # CNA projection: 1 SAMPLE 2 CHROM 3 POS 4 END 5 WAKHAN_TYPE 6 tier 7 Candidate_genes
        FILENAME == cnaf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            ss = san(\$1); sets[ss] = 1
            c = norm(\$2)
            if (!ok(c)) { offkar[ss]++; next }
            if (!isnum(\$3) || !isnum(\$4)) { badbp[ss]++; next }
            ty = toupper(\$5)
            if (ty == "GAIN" || ty == "AMP") col = "gain"
            else if (ty == "LOSS" || ty == "DEL") col = "loss"
            else { badcna[ss]++; next }
            s = \$3 + 0; e = \$4 + 0; if (e < s) { t = s; s = e; e = t }
            line = sprintf("hs%s %d %d color=%s", c, s, e, col)
            print line > "cohort.cna.txt"; print line > (ss ".cna.txt")
            ntiles["cohort"]++; ntiles[ss]++
            if (tierof(\$6) == 1 && \$7 != "" && \$7 != "NA") {
                mid = int((s + e) / 2)
                n = split(\$7, g, ";")
                for (i = 1; i <= n; i++) { gsub(/^ +| +\$/, "", g[i]); emit_labels(c, mid, g[i]) }
            }
            next
        }

        # SNV projection: 1 SAMPLE 2 CHROM 3 POS 4 SYMBOL 5 tier
        FILENAME == snvf {
            if (FNR == 1) next
            if (\$1 == "") { nosample++; next }
            ss = san(\$1); sets[ss] = 1
            if (tierof(\$5) != 1) next
            c = norm(\$2)
            if (!ok(c)) { offkar[ss]++; next }
            if (!isnum(\$3)) { badbp[ss]++; next }
            emit_labels(c, \$3 + 0, \$4)
            next
        }

        END {
            OFS = "\\t"
            print "set", "links", "tiles", "labels", "skipped_ins", "dup_mates", "off_karyotype", "bad_breakpoint", "unknown_cna_type" > "counts.tsv"
            for (s in sets) {
                tI += skip_ins[s]; tD += dup[s]; tO += offkar[s]; tB += badbp[s]; tC += badcna[s]
            }
            print "cohort", nlinks["cohort"] + 0, ntiles["cohort"] + 0, nlab["cohort"] + 0, tI + 0, tD + 0, tO + 0, tB + 0, tC + 0 > "counts.tsv"
            for (s in sets)
                print s, nlinks[s] + 0, ntiles[s] + 0, nlab[s] + 0, skip_ins[s] + 0, dup[s] + 0, offkar[s] + 0, badbp[s] + 0, badcna[s] + 0 > "counts.tsv"
            printf("INFO: cohort links=%d tiles=%d labels=%d; skipped INS=%d, duplicate mates=%d, off-karyotype=%d, bad breakpoint=%d, unknown CNA type=%d\\n",
                   nlinks["cohort"] + 0, ntiles["cohort"] + 0, nlab["cohort"] + 0, tI + 0, tD + 0, tO + 0, tB + 0, tC + 0) > "/dev/stderr"
            if (nosample > 0) printf("WARN: %d rows with an empty SAMPLE were skipped\\n", nosample) > "/dev/stderr"
        }' rec.proj.tsv sv.proj.tsv cna.proj.tsv snv.proj.tsv

    rm -f sv.proj.tsv cna.proj.tsv snv.proj.tsv rec.proj.tsv samples.txt samples.raw.txt
    """

    stub:
    """
    for s in cohort stub; do : > "\${s}.links.txt"; : > "\${s}.cna.txt"; : > "\${s}.labels.txt"; done
    printf 'set\\tlinks\\ttiles\\tlabels\\tskipped_ins\\tdup_mates\\toff_karyotype\\tbad_breakpoint\\tunknown_cna_type\\ncohort\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\nstub\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\n' > counts.tsv
    """
}
