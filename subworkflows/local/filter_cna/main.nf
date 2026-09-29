// subworkflows/local/filter_cna/main.nf
//
// Filter AnnotSV-annotated Wakhan CNVs by relevance using vembrane.
// Two sequential passes.
//   Pass 1: drop ACMG benign (∉ {1, 2}). NA/absent retained via "0" default.
//   Pass 2: (uncertain/pathogenic ACMG) AND
//           (any ST19 BED overlap OR P_gain_source OR P_loss_source).
//

include { VEMBRANE_FILTER as VEMBRANE_FILTER_PASS1 } from '../../../modules/nf-core/vembrane/filter/main'
include { VEMBRANE_FILTER as VEMBRANE_FILTER_PASS2 } from '../../../modules/nf-core/vembrane/filter/main'
include { VEMBRANE_TABLE  as VEMBRANE_TABLE_PASS1  } from '../../../modules/nf-core/vembrane/table/main'
include { VEMBRANE_TABLE  as VEMBRANE_TABLE_PASS2  } from '../../../modules/nf-core/vembrane/table/main'
include { BCFTOOLS_NORM   as BCFTOOLS_NORM_CNA     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT   as BCFTOOLS_SORT_PASS1   } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT   as BCFTOOLS_SORT_PASS2   } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW   as BCFTOOLS_VIEW_PASS1   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW   as BCFTOOLS_VIEW_PASS2   } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX  as BCFTOOLS_INDEX_PASS1  } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX  as BCFTOOLS_INDEX_PASS2  } from '../../../modules/nf-core/bcftools/index/main'
include { CREATE_DATA_DICT as CREATE_DATA_DICT_CNA } from '../../../modules/local/create_data_dict/main'

workflow FILTER_CNA {

    take:
    annotsv_vcf  // channel: [ val(meta), path(vcf.gz) ] — AnnotSV-annotated, FORMAT-reheadered Wakhan VCF
    fasta        // channel: [ val(meta2), path(fasta) ] — reference FASTA (required by bcftools/norm signature)

    main:

    // InfoTuple-safe pattern: AnnotSV fields declared Number=. arrive in vembrane
    // as InfoTuple. The chain str(x or default)...split(",")[0].strip().split(".")[0]
    // normalizes InfoTuple, plain strings, and None to a bare scalar before
    // membership comparison. NA/absent ACMG is retained (default → "0", ∉ {"1","2"}).
    def pass1_expr = """'str(INFO.get("ACMG_class") or "0").replace("(","").replace(")","").replace(chr(39),"").split(",")[0].strip().split(".")[0] not in ("1", "2")'"""

    // Pass 2 ANDed clauses:
    //   A. ACMG_class ∈ {3, 4, 5, NA, ""} — pathogenic / likely-pathogenic /
    //      uncertain. NA/absent is included so uncertain calls survive.
    //   B. Any of:
    //        - any ST19 recurrent CNV BED overlap (AMP/DEL/GAIN/LOSS) via
    //          literal-key access on name_ST19_Alterations.CNV.Recurrent_*
    //        - P_gain_source set (knotAnnotSV pathogenic gain evidence)
    //        - P_loss_source set (knotAnnotSV pathogenic loss evidence)
    //        - RE_gene set (regulatory-element gene overlap)
    //   C. CN magnitude AND confidence (Wakhan FORMAT fields):
    //        - TCN <= 1 (loss) OR TCN >= 3 (gain) — at least 1-copy deviation
    //          from diploid. Drops diploid (TCN==2) and near-diploid segments
    //          that AnnotSV may have annotated through a recurrent region.
    //        - max(CNQ1, CNQ2) >= 0.7 — at least one haplotype is confidently
    //          called. Filters out low-confidence segments where Wakhan's
    //          per-haplotype CN calls are uncertain.
    //      None-safe guards: ("TCN" in FORMAT and FORMAT["TCN"][SAMPLE] is not None)
    //      before the comparison; CNQ defaults to 0 when missing/None so
    //      max(...) >= 0.7 naturally fails.
    //   D. Haplotype imbalance filter (Wakhan CN1 / CN2 FORMAT fields):
    //        - CN1 != CN2 (allelic imbalance — the hallmark of tumor CNAs
    //          that need per-haplotype calling; drops symmetric states like
    //          balanced diploid CN1=CN2=1, tetraploid CN1=CN2=2 artefacts).
    //        - OR (CN1 == 0 AND CN2 == 0) — RESCUE regions where both
    //          haplotypes report 0. In Wakhan output this typically means
    //          the segment could not be reliably phased (rather than a true
    //          biallelic deletion); we keep those rows so downstream review
    //          can decide, rather than silently dropping them.
    //      ST19 rescue: bypass the imbalance requirement for T-ALL recurrent
    //      hits. Some ST19 loci (e.g. symmetric CDKN2A/B homdel, balanced
    //      MTAP co-deletion) present as CN1==CN2 (both alleles lost equally)
    //      yet are the highest-clinical-priority events in T-ALL. Same
    //      rescue rationale as the SV pass2 filter's ST17 rescue.
    //      None-safe guards: fall back to 0 when the FORMAT slot is missing.
    //
    // Literal-key access for ST19 fields: vembrane 2.x (cyvcf2 backend) pre-
    // collects literal INFO names at parse time. The dynamic INFO.keys()
    // fallback path triggers the broken lowercase-`info` attribute access.
    //
    // InfoTuple-safe membership pattern:
    //   str(x or "").replace("(","").replace(")","").replace(chr(39),"")
    //              .strip(", .") not in ("", "NA")

    // ST19 T-ALL recurrent CNV hit predicate — factored out for reuse across
    // clause B (existing evidence-OR) and the new clause D (ST19 rescue for
    // the CN1/CN2 balance filter). Union across the 4 SVTYPE-partitioned
    // BEDs (AMP / DEL / GAIN / LOSS).
    def st19_hit = [
        'str(INFO.get("name_ST19_Alterations.CNV.Recurrent_AMP") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST19_Alterations.CNV.Recurrent_DEL") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST19_Alterations.CNV.Recurrent_GAIN") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST19_Alterations.CNV.Recurrent_LOSS") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")'
    ].join(' or ')

    def pass2_expr = """'(str(INFO.get("ACMG_class") or "").replace("(","").replace(")","").replace(chr(39),"").split(",")[0].strip().split(".")[0] in ("3", "4", "5", "NA", "")) and ((""" + st19_hit + """) or str(INFO.get("P_gain_source") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA") or str(INFO.get("P_loss_source") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA") or str(INFO.get("RE_gene") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")) and ("TCN" in FORMAT and FORMAT["TCN"][SAMPLES[0]] is not None and (FORMAT["TCN"][SAMPLES[0]] <= 1 or FORMAT["TCN"][SAMPLES[0]] >= 3)) and (max((FORMAT["CNQ1"][SAMPLES[0]] if ("CNQ1" in FORMAT and FORMAT["CNQ1"][SAMPLES[0]] is not None) else 0), (FORMAT["CNQ2"][SAMPLES[0]] if ("CNQ2" in FORMAT and FORMAT["CNQ2"][SAMPLES[0]] is not None) else 0)) >= 0.7) and (((""" + st19_hit + """)) or (((FORMAT["CN1"][SAMPLES[0]] if ("CN1" in FORMAT and FORMAT["CN1"][SAMPLES[0]] is not None) else 0) != (FORMAT["CN2"][SAMPLES[0]] if ("CN2" in FORMAT and FORMAT["CN2"][SAMPLES[0]] is not None) else 0)) or ((FORMAT["CN1"][SAMPLES[0]] if ("CN1" in FORMAT and FORMAT["CN1"][SAMPLES[0]] is not None) else 0) == 0 and (FORMAT["CN2"][SAMPLES[0]] if ("CN2" in FORMAT and FORMAT["CN2"][SAMPLES[0]] is not None) else 0) == 0)))'"""

    // Output column list for VEMBRANE_TABLE_PASS1 and VEMBRANE_TABLE_PASS2.
    // Column order MUST stay 1:1 with the --header strings in conf/modules.config. vembrane positionally maps the i-th expression to the i-th header name.
    // BPS (Wakhan Severus-ID linkage) is now safely extractable: the
    // reheader_wakhan_fmt module rewrites Wakhan's invalid `Number=0,Type=String`
    // declaration to `Number=.,Type=String`, and the chr(44).join InfoTuple
    // pattern below handles both single-Severus-ID and multi-Severus-ID cells.
    //
    // InfoTuple-safe pattern applied to every Number=. AnnotSV field. See filter_sv/main.nf for the full rationale —
    // the chr(44).join unwraps scalar Number=. single-element tuples ("('-65',)" -> "-65"), rejoins multi-value fields with commas
    //
    // Scalar-only fields stay as plain INFO.get(K, ""): SVTYPE, CALLER, ACMG_class + AnnotSV_ranking_score, HET.
    def table_cols = [
        'CHROM', 'POS',
        // CytoBand after POS per user request. AnnotSV Number=. String;
        // qmd prepends `chr<CHROM>` at display for full-locus labels.
        'chr(44).join(str(x) for x in (INFO.get("CytoBand") or []))',
        'ID', 'REF', 'ALT',
        'INFO.get("SVTYPE", "")',
        'chr(44).join(str(x) for x in (INFO.get("SVLEN") or []))',
        'chr(44).join(str(x) for x in (INFO.get("END") or []))',
        'INFO.get("CALLER", "")',
        // HET is a Flag field in Wakhan: ##INFO=<ID=HET,Number=0,Type=Flag>.
        // cyvcf2 returns Python bool (True/False) for Flags, and vembrane's
        // cyvcf2 backend then calls .split(",") on the value during table
        // serialization, which blows up with AttributeError on bools.
        // The safe pattern for Flag fields is membership: presence in INFO
        // is the flag's semantic value. Emitted as "1"/"0" so downstream
        // consumers get a clean binary indicator column.
        '"1" if "HET" in INFO else "0"',
        'chr(44).join(str(x) for x in (INFO.get("Annotation_mode") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Gene_name") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Gene_count") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Closest_left") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Closest_right") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Location") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Overlapped_CDS_percent") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Frameshift") or []))',
        'chr(44).join(str(x) for x in (INFO.get("RE_gene") or []))',
        'INFO.get("ACMG_class", "")',
        'chr(44).join(str(x) for x in (INFO.get("ACMG") or []))',
        'INFO.get("AnnotSV_ranking_score", "")',
        'chr(44).join(str(x) for x in (INFO.get("AnnotSV_ranking_criteria") or []))',
        'chr(44).join(str(x) for x in (INFO.get("P_gain_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("P_gain_phen") or []))',
        'chr(44).join(str(x) for x in (INFO.get("P_loss_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("P_loss_phen") or []))',
        'chr(44).join(str(x) for x in (INFO.get("P_ins_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("B_gain_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("B_loss_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST19_Alterations.CNV.Recurrent_AMP") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST19_Alterations.CNV.Recurrent_DEL") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST19_Alterations.CNV.Recurrent_GAIN") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST19_Alterations.CNV.Recurrent_LOSS") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST19_Alterations.CNV.Recurrent_AMP") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST19_Alterations.CNV.Recurrent_DEL") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST19_Alterations.CNV.Recurrent_GAIN") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST19_Alterations.CNV.Recurrent_LOSS") or []))',
        'INFO.get("WAKHAN_ID", "")',
        // BPS: comma-joined Severus IDs supporting this segment's breakpoints
        // (e.g. "severus_DUP37" or "severus_BND18053,severus_BND18176").
        // Empty when Wakhan called the segment on depth only. Used downstream
        // in sv_report_index.qmd for SV-refined CNA-evidence tiering.
        'chr(44).join(str(x) for x in (INFO.get("BPS") or []))',
        'FORMAT["GT"][SAMPLE] if "GT" in FORMAT else ""',
        'FORMAT["TCN"][SAMPLE] if "TCN" in FORMAT else ""',
        'FORMAT["CN1"][SAMPLE] if "CN1" in FORMAT else ""',
        'FORMAT["CN2"][SAMPLE] if "CN2" in FORMAT else ""',
        'FORMAT["CNQ1"][SAMPLE] if "CNQ1" in FORMAT else ""',
        'FORMAT["CNQ2"][SAMPLE] if "CNQ2" in FORMAT else ""',
        'FORMAT["COV1"][SAMPLE] if "COV1" in FORMAT else ""',
        'FORMAT["COV2"][SAMPLE] if "COV2" in FORMAT else ""'
    ]
    def table_expr = table_cols.join(', ')

    // ── Split multi-allelic records ───────────────────────────────────────
    // `-m -any` splits every multi-allelic without left-alignment. nf-core module requires a fasta channel
    BCFTOOLS_NORM_CNA(
        annotsv_vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        fasta
    )

    // ── Pass 1 ───────────────────────────────────────────────────────────────
    VEMBRANE_FILTER_PASS1( BCFTOOLS_NORM_CNA.out.vcf, pass1_expr )
    BCFTOOLS_SORT_PASS1(   VEMBRANE_FILTER_PASS1.out.vcf )

    BCFTOOLS_SORT_PASS1.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_view_pass1 }

    BCFTOOLS_VIEW_PASS1(  ch_view_pass1, [], [], [] )
    BCFTOOLS_INDEX_PASS1( BCFTOOLS_VIEW_PASS1.out.vcf )
    VEMBRANE_TABLE_PASS1( BCFTOOLS_VIEW_PASS1.out.vcf, table_expr )

    // ── Pass 2 ───────────────────────────────────────────────────────────────
    VEMBRANE_FILTER_PASS2( BCFTOOLS_VIEW_PASS1.out.vcf, pass2_expr )
    BCFTOOLS_SORT_PASS2(   VEMBRANE_FILTER_PASS2.out.vcf )

    BCFTOOLS_SORT_PASS2.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_view_pass2 }

    BCFTOOLS_VIEW_PASS2(  ch_view_pass2, [], [], [] )
    BCFTOOLS_INDEX_PASS2( BCFTOOLS_VIEW_PASS2.out.vcf )
    VEMBRANE_TABLE_PASS2( BCFTOOLS_VIEW_PASS2.out.vcf, table_expr )

    // Data dictionary — same cohort-scoped pattern as FILTER_SV: `.first()`
    // + canonical meta rewrite so the parser runs exactly once (the pass2
    // CNA header schema is invariant across samples).
    CREATE_DATA_DICT_CNA(
        BCFTOOLS_VIEW_PASS2.out.vcf
            .first()
            .map { _meta, vcf -> [ [ id: 'cohort' ], vcf ] }
    )

    emit:
    vcf_pass1 = BCFTOOLS_VIEW_PASS1.out.vcf      // [ meta, *.cna.filtered.pass1.vcf.gz ]
    tbi_pass1 = BCFTOOLS_INDEX_PASS1.out.tbi
    vcf_pass2 = BCFTOOLS_VIEW_PASS2.out.vcf      // [ meta, *.cna.filtered.pass2.vcf.gz ]
    tbi_pass2 = BCFTOOLS_INDEX_PASS2.out.tbi
    tsv_pass1 = VEMBRANE_TABLE_PASS1.out.table   // [ meta, *.cna.filtered.pass1.tsv ]
    tsv_pass2 = VEMBRANE_TABLE_PASS2.out.table   // [ meta, *.cna.filtered.pass2.tsv ]
    dict      = CREATE_DATA_DICT_CNA.out.dict    // [ meta, *.data_dict.tsv ]
    versions  = channel.empty()
}
