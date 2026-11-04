# genomic_viz fixtures

Trimmed, synthetic two-sample (colo829, coloSAV) versions of the report's pass-1 cohort
tables, used by the nf-tests of `COHORT_TO_CIRCOS`, `COHORT_TO_MAF` and `MAFTOOLS_ONCOPLOT`.
Only the columns the converters read are present (plus `Gene_name_display` in the CNA table, kept for realism); the real tables have ~100. Chromosome
naming mirrors the real tables: SV/CNA `CHROM` without `chr`; `breakpoint`, `ALT` and SNV
`CHROM` with it.

## Rows and what each exercises

`cohort_sv_pass2.tsv`
- `severus_BND_ST2671_1/_2` (colo829, T1): interchromosomal BND pair, `_1` first; left gene
  overlap (BCL11B), right gene nearest (TLX3) so only BCL11B is labelled. In the map with
  n_any = 2 → thickness 4p.
- `severus_BND3978_2/_1` (colo829, T1): BND pair listed `_2` FIRST, as the real table often
  is (it is sorted by CHROM as text). The converters keep the first mate seen, so the link
  runs chr2 → chr21 and in the MAF the left gene AFF3 sits at the kept mate's locus (chr2:99727265, REF T) and the right gene RUNX1 at the partner locus (chr21:34816131, REF N); both genes overlap (AFF3,
  RUNX1). Both mates are in the map → 4p.
- `severus_DEL3918` (colo829, T2): intrachromosomal DEL; `Candidate_genes` = PCBP3 → the only
  non-BND MAF row (`SV_DEL` / `DEL`); absent from the map → 3p.
- `severus_INS13400` (colo829, T3): INS → skipped by COHORT_TO_CIRCOS by type; COHORT_TO_MAF would
  map it to `SV_INS` but writes no row only because `Candidate_genes` is `NA`.
- `severus_BND3977_1` (coloSAV, T1): unpaired `_1` (mate filtered out upstream); genes RUNX1,
  TLX3; n_any = 1 → 2p.
- `severus_INV4029` (coloSAV, T3): INV link; empty `Candidate_genes` (the real table leaves it
  blank) → no MAF row.
- `severus_BND_X1_1` (coloSAV, T3): partner on chrM → off-karyotype, dropped; `gene_name_right`
  is `unknown_contig`, as the report writes it.

`cohort_cna_pass2.tsv`: two GAIN segments (colo829; the T1 one labels OLIG2 and IFNAR1 at the
segment midpoint), one LOSS (coloSAV), one `cnloh` (unknown type → skipped and counted).

`cohort_snv.tsv`: `stop_gained` (Nonsense_Mutation / SNP); `NA` symbol (skipped);
`missense_variant&splice_region_variant` (first term wins; T1 → label NOTCH1); frameshift
deletion (Frame_Shift_Del / DEL); `weird_term` (unmapped → `unmapped_consequences.tsv`).

`cohort_sv_recurrent.map.tsv`: a Jasmine cluster map as `JASMINE_CLUSTER_MAP` writes it
(cluster prefix = index of the representative sample: `0_` colo829, `1_` coloSAV; `supp_vec`
`11` / `01`; both `BND3978` mates are listed, the ST2671 `_2` mate is not because it is dropped
as a duplicate before any lookup). Only `source_id` (= `<SAMPLE>_<ID>`) and `n_any` are read.

`candidate_genes.txt`: 12-gene panel. TMPRSS15 is deliberately NOT on it; ETV6 and PCBP3 are
(ETV6 never survives to the MAF, PCBP3 only through the DEL).

## Expected converter outputs (asserted by the nf-tests)

`COHORT_TO_CIRCOS` — links / tiles / labels: cohort 5 / 3 / 7, colo829 3 / 2 / 5,
coloSAV 2 / 1 / 3; skipped INS 1, duplicate mates 2, off-karyotype 1, unknown CNA type 1.
Labels: BCL11B, AFF3, RUNX1 (the second RUNX1 at 34.80 Mb is deduplicated into the 34 Mb
bin), TLX3, OLIG2, IFNAR1, NOTCH1.

`COHORT_TO_MAF` — 9 MAF rows (3 SNV + 5 `SV_BND` + 1 `SV_DEL`; the LEF1 deletion in MAF form chr1:101-101 T/-), 6 CN rows, 2 clinical rows,
1 unmapped term; gene list top 5 = BACH1, RUNX1, SCAF4, AFF3, BCL11B.

`mini.*` — direct `MAFTOOLS_ONCOPLOT` inputs: what `COHORT_TO_MAF` produced from the
original fixtures, before the PCBP3 DEL row was added and before the BND3978 mates were
swapped (so `mini.maf` carries the `_1` mate's chr21:34816131 / REF A coordinates and lists
SVs before SNVs). It is deliberately NOT regenerated: the onco matrix stays 10 genes ×
2 samples. Do not diff it against converter output.

## CIRCOS module inputs

`circos_three.{links,cna,labels}.txt` (two links: tier1 interchromosomal, tier2 intrachromosomal two-point form; a gain and a loss tile; two labels), `circos_single.cna.txt` (one gain tile) and `empty.txt` (zero bytes, reused for any empty track) are the direct inputs of `modules/local/circos/tests/main.nf.test`.

## COHORT_TO_CIRCOS edge cases

`c2c_edge.sv.tsv`: a sample name needing sanitisation (`P 1/a` → set `P_1_a`); a Tier-1 BND with no overlapping gene on either side (label falls back to the first `Gene_name_display` token, GENEA); a DEL with `END = NA` (counted as `bad_breakpoint`, not drawn). `c2c_edge.cna.tsv`: header only. `c2c_bad.sv.tsv`: three columns, exercises the fail-fast column check. `c2c_cohort_name.sv.tsv`: a sample literally named `cohort`, which must be rejected because it collides with the cohort track set.

## COHORT_TO_MAF edge cases

`c2m_edge.sv.tsv`: a lower-case `tra` BND pair (TRA → BND; the `_2` mate is a duplicate), an intragenic BND (same gene both sides → one row), a lower-case `dup` with a repeated gene in `Candidate_genes` (two rows, de-duplicated within the record), a row with an empty SAMPLE, a DEL with `POS = NA` (bad coordinate), and GENEQ deleted in both samples. `c2m_edge.snv.tsv`: a frameshift insertion (`Frame_Shift_Ins`, MAF form `-`/`TT` at 700-701), a +3 `protein_altering_variant` (`In_Frame_Ins`), a DNP, `splice_donor_5th_base_variant&intron_variant` (→ `Splice_Region`, not counted as altered), an `intron_variant`, an empty consequence (→ `(empty)` in the unmapped file), a `POS = NA` row, GENEK mutated in both samples, and a deletion with no shared anchor base (`ATG > C`: Start = POS, no shift). `c2m_edge.cna.tsv`: GENEP gained and lost in S1 (conflict counted, first state kept), a lower-case `amp`. `c2m_edge.panel.txt`: GENEB and GENEQ, so the gene universe with `top = 2` tests the count-desc → panel-first → name order (GENEQ before GENEK) and the cap. The SV table has no `sex`/`etp` columns, so the clinical table is header-only.

## MAFTOOLS_ONCOPLOT empty inputs

`onco_empty.maf` (header only), `onco_empty.cn.tsv`, `onco_empty.clinical.tsv` (headers only) and `empty.txt` (as the gene list) drive the `NO_DATA` paths of `modules/local/maftools/oncoplot/tests/main.nf.test`; `onco_noncoding.maf` (two rows, `Silent` and `Intron`, no oncoprint-eligible class) drives the non-coding-only `NO_DATA` path, which `read.maf` would otherwise turn into an error. `one_gene.txt` (one symbol: oncoplot needs two rows → `NO_DATA`), `mini.clinical_numeric.tsv` (sex 1/2, etp 0/1: must stay categorical), `onco_single.maf` + `onco_single.clinical.tsv` (a one-sample cohort with metadata: clinical tracks are skipped because maftools 2.26.0 collapses a one-row annotation matrix) and `onco_seven.maf` + `onco_seven.clinical.tsv` (seven samples with seven `sex` levels: more than the palette, the task must fail naming the feature) cover the same module.

## GENOMIC_VIZ subworkflow

`gv_header_only.sv.tsv`: the `cohort_sv_pass2.tsv` header with no rows. With every optional input `[]` it drives the cohort-only test of `subworkflows/local/genomic_viz/tests/main.nf.test`: one circos set whose converter outputs arrive as bare paths rather than lists, and an oncoprint `NO_DATA` marker that must still reach `viz_files`.
