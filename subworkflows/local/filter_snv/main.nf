// subworkflows/local/filter_snv/main.nf
//
// Single-pass vembrane filter for somatic SNVs.
//
// Structure:
//   (A or B or C or D) AND gnomADg_AF-rarity-gate
//
//   A. CLNSIG pathogenic-leaning (clinical override; gene-agnostic)
//   B. SCI Tier I/II somatic clinical impact (2024 fields)
//   C. ONC oncogenic (2024 fields)
//   D. Novel-variant path: HIGH/MODERATE IMPACT + both callers' VAFs above
//      threshold + not benign + SIFT deleterious AND PolyPhen damaging
//      + on candidate gene panel
//
//   VAF gate: min() over whichever per-caller VAF fields are present, so a
//   single-caller rescued variant is judged on its supporting caller rather
//   than failing because the absent caller's field defaults to 0.0.
//
// Hard gnomAD rarity filter (top-level AND):
//   `(ANN.get("gnomADg_AF") or 0.0) < snv_gnomad_af_max` — the `or 0.0`
//   fallback treats missing / blank / None AF as rare (0.0 passes any
//   positive threshold). Applied to EVERY passing variant, including
//   ClinVar-pathogenic and Tier-I SCI hits — a common ClinVar variant
//   (frequent in the population) shouldn't dominate a somatic SNV list
//   even if its ClinVar assertion is Pathogenic.
//
// SIFT + PolyPhen functional-impact consensus (clause D only):
//   Requires `"deleterious" in ANN.SIFT` AND `"damaging" in ANN.PolyPhen`
//   for the novel-variant path — a defensible tightening for variants
//   with no clinical database evidence. Kept out of A/B/C because splice
//   / nonsense / structural variants lack SIFT/PolyPhen scores.
//
// Candidate gene panel constraint is in clause D only. If params.candidate_genes
// is null/empty, the gene predicate is dropped from clause D at Groovy-string
// assembly time so vembrane doesn't error on an undefined AUX key.

include { VEMBRANE_FILTER_WITH_AUX } from '../../../modules/local/vembrane/filter_with_aux/main'
include { VEMBRANE_TABLE as VEMBRANE_TABLE_SNV } from '../../../modules/nf-core/vembrane/table/main'
include { BCFTOOLS_SORT  as BCFTOOLS_SORT_FILTER_SNV } from '../../../modules/nf-core/bcftools/sort/main'
include { CREATE_DATA_DICT as CREATE_DATA_DICT_SNV  } from '../../../modules/local/create_data_dict/main'

workflow FILTER_SNV {
    take:
    annotated_vcf     // channel: [ val(meta), path(annotated.snv.vcf.gz) ]
    annotated_tbi     // channel: [ val(meta), path(annotated.snv.vcf.gz.tbi) ]
    candidate_genes   // value:   path(candidate_genes.txt) or []

    main:

    def gene_clause = candidate_genes ?
        'and ANN.get("SYMBOL", "") in AUX["candidate_genes"]' :
        ''

    def snv_min_vaf       = params.snv_min_vaf       ?: 0.05
    def snv_gnomad_af_max = params.snv_gnomad_af_max ?: 0.01

    def filter_expr = """'(
        (
            str(INFO.get("CLNSIG") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") in ("Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance", "Likely_pathogenic_low_penetrance", "Drug_response")
            or
            str(INFO.get("SCI") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .").startswith("Tier_I")
            or
            str(INFO.get("ONC") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") in ("Oncogenic", "Likely_oncogenic")
            or
            (
                ANN.get("IMPACT", "") in ("HIGH", "MODERATE")
                and min([v for v in (INFO.get("CLAIRSTO_VAF"), INFO.get("DEEPSOMATIC_VAF")) if v is not None] or [0.0]) >= ${snv_min_vaf}
                and str(INFO.get("CLNSIG") or "").replace("(","").replace(")","").replace(chr(39),"").strip(", .") not in ("Benign", "Likely_benign")
                and "deleterious" in str(ANN.get("SIFT") or "").lower()
                and "damaging" in str(ANN.get("PolyPhen") or "").lower()
                ${gene_clause}
            )
        )
        and (ANN.get("gnomADg_AF") or 0.0) < ${snv_gnomad_af_max}
    )'""".stripIndent().replaceAll('\\n', ' ').replaceAll(' +', ' ')

    def table_expr = '''CHROM, POS, ID, REF, ALT, INFO.get("CALLER", ""), INFO.get("CLAIRSTO_VAF", ""), INFO.get("CLAIRSTO_DP", ""), INFO.get("DEEPSOMATIC_VAF", ""), INFO.get("DEEPSOMATIC_DP", ""), ANN.get("Consequence", ""), ANN.get("IMPACT", ""), ANN.get("SYMBOL", ""), ANN.get("BIOTYPE", ""), ANN.get("Feature", ""), ANN.get("HGVSc", ""), ANN.get("HGVSp", ""), ANN.get("SIFT", ""), ANN.get("PolyPhen", ""), ANN.get("gnomADg_AF", ""), ANN.get("gnomADe_AF", ""), INFO.get("CLNSIG", ""), INFO.get("CLNDN", ""), INFO.get("CLNREVSTAT", ""), INFO.get("SCI", ""), INFO.get("SCIDN", ""), INFO.get("SCIREVSTAT", ""), INFO.get("ONC", ""), INFO.get("ONCDN", ""), INFO.get("ONCREVSTAT", ""), INFO.get("ONCCONF", ""), INFO.get("GENEINFO", ""), INFO.get("RS", ""), INFO.get("VT", ""), INFO.get("CIVIC", "")'''

    // Filter
    vembrane_input_ch = annotated_vcf.join( annotated_tbi, by: 0 )
    VEMBRANE_FILTER_WITH_AUX(
        vembrane_input_ch,
        filter_expr,
        candidate_genes ?: [],
        'candidate_genes'
    )

    // Sort + bgzip + index — matches FILTER_SV's post-vembrane pattern.
    // VEMBRANE_FILTER_WITH_AUX emits [meta, vcf] (uncompressed); sort handles
    // -Oz + --write-index=tbi in ext.args (see conf/modules.config).
    BCFTOOLS_SORT_FILTER_SNV( VEMBRANE_FILTER_WITH_AUX.out.vcf )

    // Table — nf-core vembrane/table takes [meta, vcf] (no tbi) and emits .table
    VEMBRANE_TABLE_SNV( BCFTOOLS_SORT_FILTER_SNV.out.vcf, table_expr )

    // Data dictionary — cohort-scoped: run once per pipeline run using
    // `.first()` on the per-sample channel. The filtered-SNV VCF header
    // schema is stable across samples so a single dictionary suffices.
    CREATE_DATA_DICT_SNV(
        BCFTOOLS_SORT_FILTER_SNV.out.vcf
            .first()
            .map { _meta, vcf -> [ [ id: 'cohort' ], vcf ] }
    )

    emit:
    vcf  = BCFTOOLS_SORT_FILTER_SNV.out.vcf
    tbi  = BCFTOOLS_SORT_FILTER_SNV.out.tbi
    tsv  = VEMBRANE_TABLE_SNV.out.table
    dict = CREATE_DATA_DICT_SNV.out.dict
}
