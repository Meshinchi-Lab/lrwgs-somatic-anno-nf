//
// Filter AnnotSV-annotated SVs by clinical relevance using vembrane — two sequential passes.
// Pass 1: broad quality gate (precision/consensus + drop ACMG benign + drop knotAnnotSV benign gains).
// Pass 2: prioritization gate (consensus/uncertain-ACMG) AND (ST17 overlap OR pathogenic gain OR RE_gene).
//

include { VEMBRANE_FILTER     as VEMBRANE_FILTER_PASS1      } from '../../../modules/nf-core/vembrane/filter/main'
include { VEMBRANE_FILTER     as VEMBRANE_FILTER_PASS2      } from '../../../modules/nf-core/vembrane/filter/main'
include { VEMBRANE_TABLE      as VEMBRANE_TABLE_PASS1       } from '../../../modules/nf-core/vembrane/table/main'
include { VEMBRANE_TABLE      as VEMBRANE_TABLE_PASS2       } from '../../../modules/nf-core/vembrane/table/main'

include { STRIP_FLAG_CSQ                                  } from '../../../modules/local/strip_flag_csq/main'
include { BCFTOOLS_RESCUE_BND as BCFTOOLS_RESCUE_BND_PASS1 } from '../../../modules/local/bcftools/rescue_bnd/main'
include { BCFTOOLS_RESCUE_BND as BCFTOOLS_RESCUE_BND_PASS2 } from '../../../modules/local/bcftools/rescue_bnd/main'
include { BCFTOOLS_SORT       as BCFTOOLS_SORT_PASS1        } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT       as BCFTOOLS_SORT_PASS2        } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW       as BCFTOOLS_VIEW_PASS1        } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW       as BCFTOOLS_VIEW_PASS2        } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX      as BCFTOOLS_INDEX_PASS1       } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX      as BCFTOOLS_INDEX_PASS2       } from '../../../modules/nf-core/bcftools/index/main'
include { CREATE_DATA_DICT    as CREATE_DATA_DICT_SV        } from '../../../modules/local/create_data_dict/main'

workflow FILTER_SV {

    take:
    annotsv_vcf  // channel: [ val(meta), path(vcf) ]      — AnnotSV-only VCF; used when VEP is not configured
    vep_vcf      // channel: [ val(meta), path(vcf.gz) ]   — VEP VCF (AnnotSV INFO + gnomAD/ClinVar CSQ)

    main:

    // `sv_anno_vep` is the explicit opt-in for VEP on the SV branch.
    def use_vep = params.sv_anno_vep?.toString()?.toBoolean()
    // VEP can emit flag-style CSQ (";CSQ" with no value) on breakends whose ALT
    // notation it cannot parse. vembrane then treats CSQ as Flag (bool=True) and
    // any expression referencing ANN aborts with "'bool' object has no attribute
    // 'split'". STRIP_FLAG_CSQ normalizes these records so vembrane falls back
    // to its "replace with NAs" path, which the filter/table expressions handle.
    def raw_input_ch    = use_vep ? vep_vcf : annotsv_vcf
    // STRIP_FLAG_CSQ( raw_input_ch )
    // def filter_input_ch = STRIP_FLAG_CSQ.out.vcf

    // ── Pass 1: broad quality gate ───────────────────────────────────────────
    // Three hard ANDs:
    //   1. PRECISE flag set OR CALLER in (SAVANA, MINDA). SAVANA does not emit
    //      the PRECISE flag, and MINDA is the Severus+SAVANA consensus caller —
    //      both bypass the PRECISE requirement.
    //   2. ACMG_class not benign (∉ {1, 2}). NA/absent ACMG is retained because
    //      the InfoTuple-safe default coerces missing → "0", which is also
    //      outside {"1", "2"}.
    //
    // InfoTuple-safe pattern: AnnotSV fields declared Number=. arrive in
    // vembrane as InfoTuple (a Python tuple subclass). The chain
    //   str(x or default).replace("(","").replace(")","").replace(chr(39),"")
    //                    .split(",")[0].strip()
    // normalizes InfoTuple ("('3',)"), plain strings ("3", "3.0"), and None to a scalar before categorical or `.split(".")` comparison.
    // gnomAD_SV_AF max-across-InfoTuple clause. bcftools +split-vep declares  Number=.,Type=Float so vembrane wraps the value in an InfoTuple (a tuplesubclass) for a single-overlap or (0.01, 0.05, 0.02)for multi-overlap. 
    // Iterate the tuple directly and take max as float `or []` handles both "field missing" and VEP-off cases; `default=0.0` for empty tuples.
    // Interchromosomal-BND bypass for the gnomAD_SV_AF gate.
    //
    // VEP `--custom type=overlap,overlap_cutoff=80` (default `reciprocal=0`)
    // computes overlap% as `100 × intersection_bp / annotation_length_bp`. For interchromosomal BND breakpoints VEP treats the query as `ve = vs` (length 1 — the point-position); the formula SHOULD reject overlaps
    // with any large gnomAD annotation (100 / annotation_length ≪ 80). Empirically the pipeline observes 5+ high-AF gnomAD overlaps per interchromosomal BND — a known VEP quirk
    //
    // Fix: bypass the gnomAD_SV_AF gate only for INTERCHROMOSOMAL BND records — SVLEN=0 point-position queries.
    // Intrachromosomal BNDs (which can represent local DEL/INV/DUP-like rearrangements have real SVLEN>0 and VEP's 80% cutoff applies
    //
    // Normalizes an existing "chr" prefix on CHROM (some VCFs use "14", others "chr14"); 
    // the ALT bracket form always uses "chr<X>:" for the partner coord. ST17-hit predicate — factored out here and reused across (a) pass1
    // clause 3 as a rescue for high-gnomAD-AF SVs that are empirically T-ALL recurrent, (b) pass2 clause B (clinical relevance OR), and (c) pass2 clause G (MAPQ+DV rescue).
    def st17_hit = [
        'str(INFO.get("name_ST17_Alterations.SV.All_BND") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST17_Alterations.SV.All_DEL") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST17_Alterations.SV.All_DUP") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST17_Alterations.SV.All_INS") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")',
        'str(INFO.get("name_ST17_Alterations.SV.All_INV") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA")'
    ].join(' or ')

    def acmg_class_rescue = 'str(INFO.get("ACMG_class") or "0").replace("(","").replace(")","").replace(chr(39),"").split(",")[0].strip().split(".")[0] in ("3", "4", "5")'

    def pass1_expr = """'(INFO.get("PRECISE", False) or INFO.get("CALLER", "") in ("SAVANA", "MINDA") or (""" + st17_hit + """) or (""" + acmg_class_rescue + """)) and str(INFO.get("ACMG_class") or "0").replace("(","").replace(")","").replace(chr(39),"").split(",")[0].strip().split(".")[0] not in ("1", "2") and ((INFO.get("SVTYPE", "") == "BND" and ("chr" + str(CHROM).replace("chr", "") + ":") not in str(ALT)) or (""" + st17_hit + """) or max((float(x) for x in (INFO.get("gnomAD_SV_AF") or []) if x is not None and str(x) != "."), default=0.0) < 0.10)'"""

    // Table expression: single, unified across VEP-on / VEP-off. gnomAD_SV + gnomAD_SV_AF + ClinVar_SV_CLNSIG columns are always emitted; empty
    // strings when VEP off, populated by SPLIT_VEP_TO_INFO's INFO/* tags when VEP on.
    //
    // InfoTuple-safe pattern applied to every Number=. field. `INFO.get(K)` returns a pysam InfoTuple (a Python tuple subclass) for Number=. fields
    // even when only one value is present; `str()` on it yields Python's tuple repr ("('-65',)") rather than the value. The `chr(44).join(str(x) for x in (INFO.get(K) or []))`
    // iterates the tuple and rejoins with commas —  genuinely multi-value fields it preserves them as a comma-separated string.
    //
    // AnnotSV_ranking_criteria: pysam splits Number=. String fields at every comma, so joining back with "," after the split,
    // reconstructing the original pipe-separated criteria list without data loss.
    //
    // Scalar-only fields stay as plain INFO.get(K, ""). FORMAT hVAF: SEVERUS emits hVAF as three comma-separated values (HP1,HP2,unphased). REHEADER_SAVANA_FMT now declares it Number=3,
    //
    // ALT column is written as symbolic for non-BND non-INS (already <DEL>/
    // <DUP>/<INV>), as bracket notation for BND (preserved), and normalized to "<INS>" for INS records — raw INS ALT is the full inserted sequence
    // DETAILED_TYPE is from Severus (tandem_duplication, reciprocal_inv, dup_inv_segment, foldback, BFB_foldback, Reciprocal_tra, inv_tra, Templated_ins, Intra_chr_ins)
    def table_cols = [
        'CHROM', 'POS',
        // CytoBand — AnnotSV Number=. String. Prepend `chr<CHROM>` at the
        // qmd display layer so the reader sees `chr9p21` instead of a bare
        // `p21` band label.
        'chr(44).join(str(x) for x in (INFO.get("CytoBand") or []))',
        'ID', 'REF',
        '"<INS>" if INFO.get("SVTYPE") == "INS" else str(ALT)',
        'INFO.get("SVTYPE", "")',
        'INFO.get("DETAILED_TYPE", "")',
        'chr(44).join(str(x) for x in (INFO.get("SVLEN") or []))',
        'chr(44).join(str(x) for x in (INFO.get("END") or []))',
        'INFO.get("CALLER", "")',
        'INFO.get("PRECISE", "")',
        'INFO.get("SAVANA_SVTYPE", "")',
        'chr(44).join(str(x) for x in (INFO.get("Annotation_mode") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Gene_name") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Gene_count") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Closest_left") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Closest_right") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Location") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Overlapped_CDS_percent") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Frameshift") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Dist_nearest_SS") or []))',
        'chr(44).join(str(x) for x in (INFO.get("Nearest_SS_type") or []))',
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
        'chr(44).join(str(x) for x in (INFO.get("P_ins_phen") or []))',
        'chr(44).join(str(x) for x in (INFO.get("B_gain_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("B_loss_source") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST17_Alterations.SV.All_BND") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST17_Alterations.SV.All_DEL") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST17_Alterations.SV.All_DUP") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST17_Alterations.SV.All_INS") or []))',
        'chr(44).join(str(x) for x in (INFO.get("name_ST17_Alterations.SV.All_INV") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST17_Alterations.SV.All_BND") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST17_Alterations.SV.All_DEL") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST17_Alterations.SV.All_DUP") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST17_Alterations.SV.All_INS") or []))',
        'chr(44).join(str(x) for x in (INFO.get("score_ST17_Alterations.SV.All_INV") or []))',
        'FORMAT["GT"][SAMPLE] if "GT" in FORMAT else ""',
        'FORMAT["VAF"][SAMPLE] if "VAF" in FORMAT else ""',
        'chr(44).join(str(x) for x in ((FORMAT["hVAF"][SAMPLE] if "hVAF" in FORMAT else None) or []))',
        'FORMAT["DR"][SAMPLE] if "DR" in FORMAT else ""',
        'FORMAT["DV"][SAMPLE] if "DV" in FORMAT else ""',
        'chr(44).join(str(x) for x in (INFO.get("gnomAD_SV") or []))',
        'chr(44).join(str(x) for x in (INFO.get("gnomAD_SV_AF") or []))',
        'chr(44).join(str(x) for x in (INFO.get("ClinVar_SV_CLNSIG") or []))',
        'INFO.get("TUMOUR_READ_SUPPORT", "")',
        'chr(44).join(str(x) for x in (INFO.get("TUMOUR_AF") or []))',
        'INFO.get("CLASS", "")'
    ]
    def table_expr = table_cols.join(', ')

    VEMBRANE_FILTER_PASS1( raw_input_ch, pass1_expr )

    VEMBRANE_FILTER_PASS1.out.vcf
        .join( raw_input_ch, by: 0 )
        .set { ch_rescue_pass1 }

    BCFTOOLS_RESCUE_BND_PASS1( ch_rescue_pass1 )
    BCFTOOLS_SORT_PASS1( BCFTOOLS_RESCUE_BND_PASS1.out.vcf )

    BCFTOOLS_SORT_PASS1.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_view_pass1 }

    BCFTOOLS_VIEW_PASS1( ch_view_pass1, [], [], [] )
    BCFTOOLS_INDEX_PASS1( BCFTOOLS_VIEW_PASS1.out.vcf )
    VEMBRANE_TABLE_PASS1( BCFTOOLS_VIEW_PASS1.out.vcf, table_expr )

    // ── Pass 2: prioritization gate ──────────────────────────────────────────
    // Applied to pass 1 output. Pass 2 is always a strict subset of pass 1.
    //
    // Structure (all hard ANDs):
    //   A. Consensus / uncertain-ACMG   — CALLER == MINDA OR ACMG ∈ {3,4,5,NA,""}
    //   B. Clinical relevance           — ST17 hit (5 SVTYPE cols) OR P_gain
    //                                     OR RE_gene OR ClinVar pathogenic
    //   C. Rare in gnomAD               — max(gnomAD_SV_AF) < 1%
    //   D. Not inside a VNTR            — INSIDE_VNTR != "TRUE" (Severus flag,
    //                                     raw VCF value is uppercase "TRUE")
    //   E. Not in ENCODE blacklist      — both ENCODE_blacklist_left AND
    //                                     ENCODE_blacklist_right empty ("." /
    //                                     "NA" / missing)
    //   F. Phased to a single haplotype — INFO/HP ∈ (1, 2). Somatic sanity
    //                                     filter: HP=0 or missing = unphased,
    //                                     less confident single-haplotype origin.
    //   G. Read-support quality gate    — (MAPQ >= 20 AND DV >= 5) OR ST17 hit.
    //                                     ST17-hit records rescue borderline
    //                                     MAPQ / DV — a T-ALL recurrent SV in
    //                                     ST17 is high-value even at lower
    //                                     read depth (SVTYPE-specific column
    //                                     match; not a bare presence check).
    //
    // Each name_ST17_<stem> field is accessed by literal key rather than via
    // INFO.keys() iteration. vembrane 2.x (cyvcf2 backend) statically pre-
    // collects literal INFO field names at parse time and wires up typed
    // accessors through cyvcf2's `Variant.INFO` (uppercase). The dynamic
    // INFO.keys() fallback path in 2.4.0 tries `Variant.info` (lowercase) and
    // aborts on the first record with
    //   "'cyvcf2.cyvcf2.Variant' object has no attribute 'info'".
    // Literal-key access avoids that broken fallback entirely.
    //
    // InfoTuple-safe membership pattern (P_gain_source / RE_gene / ST17 /
    // ENCODE_blacklist_* cols):
    //   str(x or "").replace("(","").replace(")","").replace(chr(39),"")
    //              .strip(", .") not in ("", "NA")
    // → True iff the field has real annotation content (handles InfoTuple,
    //   plain string, None, and placeholders ".", "NA", "()" uniformly).
    // For blacklist columns the sense is inverted: we want "in" (empty) rather
    // than "not in" (populated) — a populated blacklist column means the SV
    // sits inside a dark region and should be dropped.
    //
    // `st17_hit` is defined above pass1_expr and reused here in pass2
    // clauses B and G (unchanged behavior — same predicate, single source
    // of truth).

    def pass2_clauses = [
        // A. Consensus / uncertain-ACMG
        '(INFO.get("CALLER", "") == "MINDA" or str(INFO.get("ACMG_class") or "").replace("(","").replace(")","").replace(chr(39),"").split(",")[0].strip().split(".")[0] in ("3", "4", "5", "NA", ""))',
        // B. Clinical relevance
        '(' + st17_hit + ' or str(INFO.get("P_gain_source") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA") or str(INFO.get("RE_gene") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("", "NA") or "athogenic" in chr(44).join(str(x) for x in (INFO.get("ClinVar_SV_CLNSIG") or [])))',
        // C. Rare in gnomAD (<1%) — bypassed for (a) interchromosomal BNDs
        //    where VEP's `--custom type=overlap,overlap_cutoff=80` formula
        //    is broken on SVLEN=0 point-position queries, or (b) any SV
        //    hitting the ST17 T-ALL recurrence BED, where empirical T-ALL
        //    recurrence overrides population-AF concerns (the 9q34 ABL1/
        //    NUP214/SET locus and BCL11B/TLX3 regions overlap common
        //    polymorphic gnomAD CNVs, but the somatic T-ALL event is
        //    unrelated to the germline polymorphism). Same pattern as
        //    pass1 clause 3 — both passes use identical rescue logic.
        '((INFO.get("SVTYPE", "") == "BND" and ("chr" + str(CHROM).replace("chr", "") + ":") not in str(ALT)) or (' + st17_hit + ') or max((float(x) for x in (INFO.get("gnomAD_SV_AF") or []) if x is not None and str(x) != "."), default=0.0) < 0.01)',
        // D. Not inside VNTR (case-insensitive; Severus emits raw "TRUE")
        'str(INFO.get("INSIDE_VNTR") or "").strip().lower() != "true"',
        // E. Neither ENCODE blacklist column populated
        'str(INFO.get("ENCODE_blacklist_left") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") in ("", "NA")',
        'str(INFO.get("ENCODE_blacklist_right") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") in ("", "NA")',
        // F. Haplotype-supported somatic signal via hVAF, with ST17 rescue.
        //    Severus hVAF header: (H0, H1, H2) where H0 is unphased, H1 is
        //    haplotype 1, H2 is haplotype 2 (see Severus FORMAT declaration
        //    ##FORMAT=<ID=hVAF,Number=3,Type=Float,Description="Haplotype
        //    specific variant Allele frequency (H0,H1,H2)">).
        //    Rule (inclusive OR): at least one of hVAF[H1] or hVAF[H2] > 0.
        //    Equivalently: NOT (H1 == 0 AND H2 == 0). Records with all-zero
        //    phased haplotype support are dropped as caller anomalies.
        //    Slice `[1:3]` deliberately excludes the unphased H0 slot — a
        //    somatic SV should be supported by at least one *phased*
        //    haplotype, not just unphased reads.
        //    The INFO/HP tag is unreliable on real Severus output — missing
        //    (`.`) on the majority of records and numerically inconsistent
        //    with hVAF slots when set — so we rely on FORMAT/hVAF directly.
        //    ST17 rescue retained: a T-ALL recurrent hit with borderline
        //    haplotype support keeps clinical priority.
        '(any((v or 0.0) > 0.0 for v in ((FORMAT["hVAF"][SAMPLES[0]] if "hVAF" in FORMAT else None) or (0.0, 0.0, 0.0))[1:3]) or ' + st17_hit + ')',
        // G. MAPQ+DV quality gate OR ST17 rescue. vembrane FILTER context
        //    exposes samples as SAMPLES[0] (indexed into the sample list),
        //    NOT SAMPLE — SAMPLE is a per-column placeholder specific to
        //    `vembrane table`. Same idiom as filter_cna's TCN/CNQ1 clauses.
        '(((INFO.get("MAPQ") or 0) >= 20 and ("DV" in FORMAT and FORMAT["DV"][SAMPLES[0]] is not None and FORMAT["DV"][SAMPLES[0]] >= 5)) or ' + st17_hit + ')'
    ]
    def pass2_expr = "'" + pass2_clauses.join(' and ') + "'"

    VEMBRANE_FILTER_PASS2( BCFTOOLS_VIEW_PASS1.out.vcf, pass2_expr )

    // BND rescue for pass 2 searches the pass 1 bgzipped output for orphaned mates.
    VEMBRANE_FILTER_PASS2.out.vcf
        .join( BCFTOOLS_VIEW_PASS1.out.vcf, by: 0 )
        .set { ch_rescue_pass2 }

    BCFTOOLS_RESCUE_BND_PASS2( ch_rescue_pass2 )
    BCFTOOLS_SORT_PASS2( VEMBRANE_FILTER_PASS2.out.vcf )

    BCFTOOLS_SORT_PASS2.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_view_pass2 }

    BCFTOOLS_VIEW_PASS2( ch_view_pass2, [], [], [] )
    BCFTOOLS_INDEX_PASS2( BCFTOOLS_VIEW_PASS2.out.vcf )
    VEMBRANE_TABLE_PASS2( BCFTOOLS_VIEW_PASS2.out.vcf, table_expr )

    // Data dictionary — parse the pass2 VCF's INFO/FORMAT/FILTER header lines
    // into a TSV keyed by field_name so SV_REPORT_INDEX can render column
    // definitions above the SV DT table. Consumed via CREATE_DATA_DICT_SV.out.dict.
    //
    // Cohort-scoped: the pass2 VCF header schema is INVARIANT across samples
    // (same source callers → same AnnotSV/VEP pipeline → same declarations),
    // so running the parser once per sample would emit N identical TSVs.
    // `.first()` collapses the per-sample channel to a single tuple, and the
    // meta rewrite to `[id: 'cohort']` yields a stable filename that isn't
    // tied to whichever sample happened to arrive first.
    CREATE_DATA_DICT_SV(
        BCFTOOLS_VIEW_PASS2.out.vcf
            .first()
            .map { _meta, vcf -> [ [ id: 'cohort' ], vcf ] }
    )

    emit:
    vcf_pass1 = BCFTOOLS_VIEW_PASS1.out.vcf           // [ meta, *.filtered.pass1.vcf.gz ]
    tbi_pass1 = BCFTOOLS_INDEX_PASS1.out.tbi           // [ meta, *.filtered.pass1.vcf.gz.tbi ]
    vcf_pass2 = BCFTOOLS_VIEW_PASS2.out.vcf           // [ meta, *.filtered.pass2.vcf.gz ]
    tbi_pass2 = BCFTOOLS_INDEX_PASS2.out.tbi           // [ meta, *.filtered.pass2.vcf.gz.tbi ]
    tsv_pass1 = VEMBRANE_TABLE_PASS1.out.table        // [ meta, *.filtered.pass1.tsv ]
    tsv_pass2 = VEMBRANE_TABLE_PASS2.out.table        // [ meta, *.filtered.pass2.tsv ]
    dict      = CREATE_DATA_DICT_SV.out.dict          // [ meta, *.data_dict.tsv ]
    versions  = channel.empty()
}
