# WGS Nanopore Annotation

A Nextflow DSL2 pipeline for merging, annotating, filtering, and prioritising **structural variants (SVs)**, **copy number alterations (CNAs)**, and **somatic single-nucleotide variants (SNVs)** from long-read Nanopore whole-genome sequencing. Developed for T-ALL cancer genomics at the Quintarelli and Locatelli laboratories at Ospedale Bambino Gesu.

## Overview

Input VCFs are produced upstream by two nf-core-style pipelines that this workflow consumes directly:

- **[KolmogorovLab/Lumos](https://github.com/KolmogorovLab/Lumos)** (run from the [Meshinchi-Lab/Lumos](https://github.com/Meshinchi-Lab/Lumos) fork) — long-read WGS SV calling (Severus), allele-specific copy number calling (Wakhan), and somatic SNV calling (DeepSomatic) from Nanopore BAMs. This cohort is **tumour-only** (no matched normal), and default parameters are based on this assumption. 
- **[epi2me-labs/wf-somatic-variation](https://github.com/epi2me-labs/wf-somatic-variation)** — Oxford Nanopore's end-to-end somatic workflow, providing ClairS-TO SNV calls.

SAVANA structural variants are called within this pipeline, unless precomputed SAVANA VCFs are supplied in the samplesheet `savana_sv` column.

The pipeline then runs three parallel variant-type paths — SV, CNA, SNV — that generates a single cohort-level Quarto report / Interactive dashboard with genomic data visualization:

![](assets/pipeline_graphic.png)

---

## Subworkflows

A parallel flow is used for all three variant types (SV, CNA, SNV): (i) merge / consensus across per-caller VCFs, (ii) annotate with knowledge bases + prediction tools, (iii) vembrane-based filtering with driver-tier prioritization, (iv) render into the cohort dashboard.

| Variant type | Merge / consensus | Annotate | Filter | Report |
|---|---|---|---|---|
| SV  | MERGE_SV + MINDA_ANNOTATIONS | ANNOTATIONS_SV  | FILTER_SV  | SV_REPORT_INDEX |
| CNA | Wakhan Only    | ANNOTATIONS_CNA | FILTER_CNA | SV_REPORT_INDEX |
| SNV | MERGE_SNV                    | ANNOTATIONS_SNV | FILTER_SNV | SV_REPORT_INDEX |

### MERGE_SV

Prepares and merges the per-caller VCFs:

1. `BCFTOOLS_REHEADER` — standardise sample headers
2. `BCFTOOLS_VIEW` / `BCFTOOLS_SORT` / `BCFTOOLS_INDEX` — filter, sort, and index each caller VCF
3. Merge with **MINDA** — ensemble merging with CALLER-aware consensus


### MINDA_ANNOTATIONS

Cross-references the MINDA ensemble VCF with the original per-caller VCFs to produce two outputs used downstream:

- **Annotated Severus VCF** — carries `CALLER`, `SAVANA_ID`, `SAVANA_SVTYPE`, `SAVANA_POS`, `SAVANA_END`, `SAVANA_CHR2`, `SAVANA_POS2`, `SAVANA_MATEID` INFO fields annotated from the MINDA overlap, plus SAVANA's own `TUMOUR_READ_SUPPORT`, `TUMOUR_AF` and `CLASS` for consensus record so the SV tables and report show SAVANA's somatic quality metrics

- **SAVANA-only VCF** — calls present in SAVANA but absent from the MINDA ensemble (sorted and indexed)

The `CALLER` field values are:
| Value | Meaning |
|-------|---------|
| `MINDA` | Call supported by both Severus and SAVANA |
| `SEVERUS` | Severus-only call |
| `SAVANA` | SAVANA-only call (in the SAVANA-only VCF) |

### ANNOTATIONS_SV

1. **Concatenate** annotated Severus VCF + SAVANA-only VCF (`BCFTOOLS_CONCAT`)
2. **Sort + index** the merged VCF
3. **AnnotSV** annotation — structural variant annotation with optional user-supplied BED files (FtIncludedInSV / SVincludedInFt / AnyOverlap overlap categories); see [User BED annotations](#user-bed-annotations) below
4. **BCFTOOLS_REHEADER_SAVANA_FMT** — restore the Severus/SAVANA FORMAT definitions AnnotSV drops, and put contig names back on the reference FASTA's, read from its `.fai` (`bcftools annotate --rename-chrs` + `bcftools reheader --fai`). AnnotSV strips `chr` from CHROM and declares an Ensembl-style header (`MT`, not `M`), which broke bcftools on mitochondrial calls and FASTA lookups downstream (Jasmine `--dup_to_ins`). A UCSC-style FASTA gives `chr1`…`chrM`, an Ensembl-style one `1`…`MT`; a CHROM with no FASTA contig stops the run with the names involved
5. **KnotAnnotSV** — interactive HTML (or `.xlsm`) report from the AnnotSV TSV output
6. **ENSEMBLVEP_VEP** (optional) — VEP annotation of the AnnotSV VCF with `--custom` gnomAD SV and ClinVar SV overlaps; requires either a local cache path (`--vep_cache`) or `include_cache = true` to auto-download. Both reference VCFs first pass through **SV_REFERENCE_FILTER_CHECK** (in `PREPARE_REFERENCES`): a FILTER-column census is published to `reference/sv_filter_check/`, and under the default `--sv_reference_filter_policy filter` VEP annotates against the PASS-only subset (gnomAD-SV v4.1 sites ships 955,369 non-PASS records out of 2,154,486; dbVar nstd102 carries no FILTER annotation and is used as-is).

### FILTER_SV

Applies a clinical relevance filter to the VEP-annotated VCF (or the AnnotSV VCF when VEP is not configured) using [vembrane](https://github.com/vembrane/vembrane).

**Hard filters** (all must pass):

| Filter | Field | Criterion |
|--------|-------|-----------|
| Precise breakpoints | `PRECISE` flag | Must be present |
| Not benign | `ACMG_class` | Not class 1 (benign) or class 2 (likely benign) |
| gnomAD population frequency | `ANN["gnomAD_SV_AF"]` (VEP CSQ) | All overlapping gnomAD SVs must have AF < 10% |

**Inclusion criteria** (at least one must be true):

| Criterion | Field |
|-----------|-------|
| Consensus call (Severus + SAVANA) | `CALLER == "MINDA"` |
| VUS or higher | `ACMG_class` 3, 4, or 5 |
| Pathogenic in AnnotSV databases | `P_gain_source`, `P_loss_source`, or `P_ins_source` non-empty |
| ClinVar pathogenic | `ANN["ClinVar_SV_CLNSIG"]` contains "Pathogenic" |
| T-ALL cohort SV overlap | `name_ST17_*` INFO field non-empty (user BED annotations) |
| COSMIC somatic SV | `Cosmic_ID` non-empty |

Additional per-record refinements: 
 * hVAF haplotype-imbalance filter
 * Tier-3 MINDA-only rescue (T3 records must have consensus support unless they overlap an ST17 BED region). 
 * When VEP is not configured the gnomAD and ClinVar checks are omitted and vembrane filters on AnnotSV INFO fields only.

Output: bgzipped, tabix-indexed filtered VCF (`*.merged.annotated.filtered.pass2.vcf.gz`) plus a vembrane-derived TSV table.

### ANNOTATIONS_CNA

Mirrors ANNOTATIONS_SV for the Wakhan CNV VCF:

1. `BCFTOOLS_REHEADER` — standardise sample header to `meta.id`
2. `BCFTOOLS_RENAME_IDS` — standardise IDs to be compatible with AnnotSV 
3. `BCFTOOLS_ADD_CALLER_TAG` — tag every record with `INFO/CALLER=WAKHAN` 
4. `BCFTOOLS_SORT` + `BCFTOOLS_INDEX`
5. `ANNOTSV_ANNOTSV` — CNV annotation with the same user-BED overlap set as SV (ST19 T-ALL recurrent-CNV BEDs live under `Users/`)
6. `BCFTOOLS_REHEADER_WAKHAN_FMT` — Fix header and dropped values from AnnotSV via `bcftools query` + `bcftools annotate` with an auditable TSV output, and put contig names back on the reference FASTA's as in ANNOTATIONS_SV. See module header for the full traceable / auditable fix.
7. `BCFTOOLS_CLEAN_ANNOTSV` — normalise `ACMG_class` and `AnnotSV_ranking_score` to typed scalars (AnnotSV emits compound `<split>,full=<full>` encodings that break vembrane)
8. `KNOTANNOTSV` — interactive HTML + Excel reports, includes all reported variants

### FILTER_CNA

Two-pass vembrane filter on the AnnotSV-annotated Wakhan VCF.

**Pass 1** — drop ACMG benign:
- Records with `ACMG_class in {1, 2}` removed; NA / absent ACMG retained 

**Pass 2** — prioritise: `(ACMG uncertain/pathogenic) AND (ST19 hit OR knotAnnotSV pathogenic evidence OR RE_gene) AND (TCN magnitude + CNQ confidence) AND (haplotype-imbalance OR ST19 rescue)`:

| Clause | Field | Criterion |
|--------|-------|-----------|
| A. ACMG bucket | `ACMG_class` | ∈ {3, 4, 5, NA} — includes uncertain |
| B. Evidence OR-gate | ST19 BED, `P_gain_source`, `P_loss_source`, `RE_gene` | at least one non-empty |
| C. TCN + confidence | `FORMAT/TCN`, `FORMAT/CNQ1`, `FORMAT/CNQ2` | TCN ≤ 1 (loss) or ≥ 3 (gain); max(CNQ1, CNQ2) ≥ 0.7 |
| D. Haplotype imbalance | `FORMAT/CN1`, `FORMAT/CN2` | CN1 ≠ CN2 OR (CN1 == 0 AND CN2 == 0) — rescues unhaplotypable regions. ST19 hit bypasses the imbalance requirement (biallelic-null CDKN2A/B homdels present as CN1==CN2==0). |

Output: bgzipped, tabix-indexed filtered VCF (`*.cna.filtered.pass2.vcf.gz`) plus a vembrane-derived TSV table including the BPS (severus IDs) column for SV↔CNA cross-reference.

### MERGE_SNV

Ensembles two long-read somatic small-variant callers:

1. `BCFTOOLS_NORM` — per-caller multi-allelic split + left-alignment + trim (ClairS-TO, DeepSomatic)
2. `BCFTOOLS_VIEW -f PASS` — PASS-filter each caller
3. `BCFTOOLS_ISEC -n+2 -c none` — strict 2-of-2 intersection between ClairS-TO and DeepSomatic 
4. `RENAME_CALLER_FIELDS` — rename `FORMAT/AF` (ClairS-TO) and `FORMAT/VAF` (DeepSomatic) to caller-tagged names to distinguish per-caller values
5. `VIEW_TO_INFO` — retain both callers' VAF/DP in the INFO tags
6. `BCFTOOLS_ANNOTATE_CROSS_CALLER` — a single row carries both callers' read support evidence
7. `BCFTOOLS_ADD_CALLER_TAG` — set `INFO/CALLER=CLAIRSTO,DEEPSOMATIC` for provenance
8. `BCFTOOLS_SORT` + `BCFTOOLS_REHEADER` — final sort; sample column renamed to `meta.id` so cohort-level joins group by sample correctly
9. `QC_MERGE_SNV` — warning-only data quality check comparing the final consensus VCF against the bcftools isec ground truth 

Output: `*.consensus.snv.vcf.gz` + tabix — the input to ANNOTATIONS_SNV.

### ANNOTATIONS_SNV

Layered annotation of the consensus SNV VCF:

1. `BCFTOOLS_ANNOTATE_CLINVAR_SNV` — copy ClinVar INFO fields (`CLNSIG`, `CLNDN`, `CLNREVSTAT`, `SCI`, `ONC`, etc.) from the normalized ClinVar VCF via `bcftools annotate -c INFO/...`. Uses SETUPCLINVAR's output (Bioconductor-normalised chr-prefixed ClinVar VCF with 2024 somatic-specific tags SCI + ONC preserved).
2. `BCFTOOLS_ANNOTATE_CIVIC_SNV` — copy `INFO/VT` (CIViC Variant Name) and `INFO/CIVIC` (40-field pipe-encoded CIViC assertion payload) from SETUPCIVIC's output. CIViC's raw INFO field is named `CSQ`; SETUPCIVIC renames it to `CIVIC` via `bcftools annotate --rename-annots` so it doesn't collide with VEP's downstream CSQ.
3. `ENSEMBLVEP_VEP_SNV` — VEP consequence + IMPACT + gnomAD frequency (`--af_gnomade --af_gnomadg`) in a fresh `INFO/CSQ` field. Runs offline with `--pick_allele_gene` (one consequence per allele+gene, canonical transcript preferred) and `--sift b --polyphen b`.
4. `BCFTOOLS_SORT_ANNOTATIONS_SNV` — final sort + tabix index.

Output: fully-annotated SNV VCF with ClinVar + CIViC + VEP annotations, feeds FILTER_SNV.

### FILTER_SNV

Single-pass vembrane filter over the fully-annotated SNV VCF.


| Clause | Criterion |
|--------|-----------|
| A. ClinVar pathogenic-leaning | `CLNSIG` matches "Pathogenic" / "Likely_pathogenic" / "Pathogenic_low_penetrance" |
| B. SCI Tier I / II somatic clinical impact | 2024 ClinVar `SCI` field not empty |
| C. ClinVar oncogenic | `ONC` contains "Oncogenic" or "Likely_oncogenic" |
| D. Novel-variant path | HIGH/MODERATE IMPACT + `CLAIRSTO_VAF` **or** `DEEPSOMATIC_VAF` ≥ `params.snv_min_vaf` (default 0.05) + not benign + SIFT `deleterious` + PolyPhen `damaging` + on `--candidateGenesFile` panel |

`ANN.gnomADg_AF < params.snv_gnomad_af_max` (default 0.01) applied to every clause — a common ClinVar-Pathogenic variant that is frequent in gnomAD is dropped as a soft PoN.

Output: bgzipped, tabix-indexed filtered VCF (`*.snv.filtered.vcf.gz`) plus a vembrane-derived TSV table with per-caller VAFs, ClinVar / CIViC fields, and VEP consequence + IMPACT + gnomAD AF.

### SV_REPORT_INDEX (cohort-level)

Collects per-sample outputs from all samples across the cohort from  `FILTER_SV`, `FILTER_CNA`, and `FILTER_SNV` modules and renders a single `index.html` interactive dashboard cohort report using [Quarto](https://quarto.org) (R/`rocker:verse` container).

Inputs gathered across all samples:
- AnnotSV annotation TSVs
- KnotAnnotSV interactive HTML reports
- VEP summary HTML reports (when VEP is configured)
- Optional sample metadata CSV (`--sample_metadata`)

The rendered `index.html` provides a cohort-level view of all per-sample IGV, KnotAnnotSV, and VEP reports in one place, suitable for cohort-level SV review and filtering.

If `run_cohort_merge_sv` the same report, the SV page gains a **Cohort recurrent SVs (Jasmine merge)** card and the tiered SV table (and the published `cohort_sv_pass2.tsv`) record five recurrence columns joined by sample and record ID: `cohort_recurrent`, `n_any`, `n_tiered`, `cluster_span`, `cluster_id`.

The current SV report is developed for T-ALL and assumes `SVinFt_DEL` and `SVinFt_DUP` are generated from BEDs are the ST17 DEL/DUP references. SV_REPORT_INDEX merges the two columns into `name_ST17_DEL`/`name_ST17_DUP` and their scores before its ST17 cleaning step, so the type matching, the BND partner check, `ST17_hit` and the driver tier all see one set of matches.

Output: `results/sv_report/index.html`

### GENOMIC_VIZ — cohort circos and oncoprint

Runs when `run_genomic_viz` is set (default true). A `CIRCOS` renders PNG + SVG for the cohort using the Circlize R package. Next, filtered TSV are converted into a **MAF** + copy-number table + gene list for **maftools**, and an oncoprint is generated with  `MAFTOOLS_ONCOPLOT` (Bioconductor 2.26.0) to create a gene × sample oncoprint with SNV classes from VEP consequence. 

`SV_*` classes and Amp/Del. Encodings: link colour = `driver_tier`, thickness = cohort `n_any`, Tier-1 gene labels (SV breakpoint genes, copy-number candidate genes, SNV genes); oncoprint rows = candidate-panel genes with a non-synonymous SNV, SV or copy-number alteration plus genes with such an alteration in ≥ 2 samples, limited by `oncoprint_top_genes` (default 40). 

Figures are included in the **Cohort figures** page: the cohort circos and oncoprint embedded with captions that state the encodings and counts, a thumbnail gallery of the per-sample circos, and every intermediate downloadable. The same files are also published under `genomic_viz/`. 

---

## Information

**Filtering strategy and driver-tier scoring**

Every variant type (SV, CNA, SNV): (1) a PASS filter at the caller / vembrane step, (2) a soft evidence-tier assignment based on empirical T-ALL recurrence + candidate-gene panel membership + curated pathogenicity, and (3) a per-tier tie-breaker `tier_score` used to rank rows inside each tier. 

Tier 1 variants are included in the interactive DISCO plot and lollipop plots. 

The `ACMG_class` used by the SV and CNA filters and tiers comes from AnnotSV, which applies the germline copy-number framework (Riggs et al. 2020). For how that relates to the somatic standard for these tumour-only calls (Mikhail et al. 2019), see [Guidelines for structural and copy-number variants](docs/acmg_sv_cnv_guidelines.md).

### Structural variants (SV) 

| Tier | Criterion | `tier_score` tie-breaker |
|---|---|---|
| **T1** | Any SVTYPE-matched `name_ST17_*` hit (empirical T-ALL recurrence BED). For DEL and DUP records these include the `SVINFT_ANNOTATE` matches. Partner chromosome is cross-validated against the ALT bracket notation within a 10 kb tolerance for interchromosomal BND types. | max integer score across the SVTYPE-matched `score_ST17_*` is an estimate for the number many T-ALL samples in the reference with this event. |
| **T2** | No T-ALL reference hit, but ≥ 1 gene from `Gene_name` / `Closest_left` / `Closest_right` overlaps the `--candidateGenesFile` panel. | count of overlapping candidate genes (breadth of impact). |
| **T3** | Neither an T-ALL reference hit nor a candidate-gene overlap, but the record has pathogenic evidence (P_gain/P_loss/P_ins source, RE_gene, ACMG 3–5, ClinVar pathogenic). | `AnnotSV_ranking_score`. |

### Copy-number alterations (CNA)

| Tier | Criterion | `tier_score` tie-breaker |
|---|---|---|
| **T1** | Any of the four SVTYPE-matched `name_ST19_Alterations.CNV.Recurrent_*` (empirical T-ALL recurrence BED) hits. GAIN records match AMP + GAIN columns; LOSS records match DEL + LOSS columns. | max integer score across the four SVTYPE-matched `score_ST19_*` |
| **T2** | Same candidate-gene rule as SV. | count of overlapping candidate genes. |
| **T3** | Passed the FILTER_CNA gate on other evidence (P_gain_source, ACMG 3–5, etc.). | `AnnotSV_ranking_score`. |

### Single-nucleotide variants (SNV)

| Tier | Criterion | `tier_score` tie-breaker |
|---|---|---|
| **T1** | ClinVar oncogenicity `ONC` field contains "oncogenic" ("Oncogenic" + "Likely_oncogenic") **AND** `SYMBOL` is in the candidate-gene `--candidateGenesFile` panel. | max(`CLAIRSTO_VAF`, `DEEPSOMATIC_VAF`) reports the higher VAF somatic call. |
| **T2** | ClinVar `CLNSIG` matches "pathogenic" (catches "Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance") **OR** `SYMBOL` is in the candidate-gene panel. | same VAF max. |
| **T3** | Has some pathogenicity evidence (SIFT deleterious / PolyPhen damaging predictions). | same VAF max. |


### Custom References 

The `--candidateGenesFile` panel and the customized BED files from the Nextflow.config are found in the `--outdir` under annotsv/setupuseranno so these can be used for further down-stream investigations and comparisons. 

---

### User BED annotations

User BED files supplied via `--FtIncludedInSV`, `--SVincludedInFt`, or `--AnyOverlap` are staged into the AnnotSV reference directory by the `ANNOTSV_SETUPUSERANNO` module before annotation runs.

**BED file requirements:**

- Each file must include a tab-delimited header comment as the first line, beginning with `#chrom`:

  ```
  #chrom<TAB>chromStart<TAB>chromEnd<TAB>name<TAB>score<TAB>strand
  ```

- Column delimiters must be **tabs** (not spaces). AnnotSV uses the column names from this line as VCF INFO field names.
- The module appends the file basename (without `.bed`) to every column name after `#chrom`, producing unique INFO fields per file (e.g. `name_ST17_Alterations.SV.All_DEL`). The filter expression matches these fields with the `name_ST17_*` pattern.
- The module validates that all `header.tsv` files are non-empty and will exit with a descriptive error if a BED file is missing its header line or uses space delimiters.

**Reference DEL/DUP matches in the SVincludedInFt method (`SVINFT_ANNOTATE`).** AnnotSV FtIncludedInSV parameter, for example, will identify a 2.1 kb BCL11B-enhancer tandem duplication falls inside the 2.6 kb ST17 entries for the same event. 

After FILTER_SV, `SVINFT_ANNOTATE` uses an R script (`bin/svinft_annotate.R`, plyranges, in the `sv_report_index` image) on the filtered prioritized SVs to compare them again against the reference DEL and DUP BEDs given by `--svinft_del_bed` and `--svinft_dup_bed` (same format as the user BEDs above). Any reference SV that falls within a sample's SV of the same type, covered by at least `--svinft_min_ref_coverage` (default 0.5). The length of coordinates overlap prevents megabase length entries from matching small SV that fall within large events. The module generates two columns, `SVinFt_DEL` and `SVinFt_DUP`, indicating if the reference SV (SV) overlaps in the sample feature (FT) coordinates. The EnsDb (`--ensembl_db`) supplies chromosome names and lengths.

The step is skipped, with a warning, when either BED or the EnsDb is unset or not found. The augmented tables are published under `filter_sv/svinft/`. 

---

## Requirements

- Nextflow 25.10.0 (pinned in `nextflow.config`)
- Docker, Singularity/Apptainer, or Conda
- Reference files (see [Reference data](#reference-data) below)

---

## Usage

### Samplesheet

The pipeline takes a CSV samplesheet. Each row is one sample with paths to aligned BAM/CRAM files and pre-called SV/CNV/SNV VCFs from the upstream data processing pipelines Lumos and wf-somatic-variation:

```csv
id,bam,severus_sv,savana_sv,wahkan_cnv,clairsto_snv,deepsomatic_snv
SAMPLE_1,/path/to/SAMPLE_1.haplotagged.bam,/path/to/SAMPLE_1.severus_somatic.vcf.gz,,/path/to/SAMPLE_1.wakhan_cna.vcf.gz,/path/to/SAMPLE_1.wf-somatic-snv.vcf.gz,/path/to/SAMPLE_1.deepsomatic.vcf.gz
```

`savana_sv` is optional. Leave it empty (as above), write `NA`, or omit the column to have the pipeline run SAVANA. To supply a precomputed call set, provide the SAVANA somatic VCF output, named `<sample>.classified.somatic.vcf` (as written by `savana to`). 

### Minimal run

```bash
nextflow run Meshinchi-Lab/lrwgs-somatic-anno-nf \
    -r main \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `data/samplesheet.csv` | Path to samplesheet CSV |
| `--outdir` | `results/` | Output directory |
| `--genome` | `GRCh38` | Genome build |
| `--fasta` | `data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | Path to reference FASTA |
| `--contigs_file` | SAVANA `contigs.chr.hg38.txt` (GitHub URL) | Chromosome contig list for SAVANA |
| `--pon_1kg` | `1000g_hg38` | SAVANA panel-of-normals identifier |
| `--sv_merge_size` | `100` | Size threshold in base pairs for SV merging across Severus and SAVANA |
| `--annotsv_annotations` | `data/reference/AnnotSV` | Path to pre-built AnnotSV database (null = auto-download) |
| `--candidateGenesFile` | `data/reference/custom/cancer_genes.txt` | Gene list for AnnotSV candidate gene annotation and used in the SV HTML Report for driver tier classification. A list of genes known to be altered in a specific cancer and/or pan-cancer associated genes |
| `--FtIncludedInSV` | `*{DEL,DUP,INS,INV}.bed` in `data/reference/custom/` | BED files for AnnotSV FtIncludedInSV user annotations |
| `--SVincludedInFt` | `[]` | BED files for AnnotSV SVincludedInFt user annotations |
| `--AnyOverlap` | `data/reference/custom/*BND.bed` | BED files for AnnotSV AnyOverlap user annotations |
| `--svinft_del_bed`, `--svinft_dup_bed` | `*_{DEL,DUP}.bed` in `data/reference/custom/` | Reference DEL and DUP BEDs for `SVINFT_ANNOTATE`; unset either one to skip the step |
| `--svinft_min_ref_coverage` | `0.5` | Fraction of a reference entry's length a DEL/DUP must cover to match it in `SVINFT_ANNOTATE`, in (0, 1] |
| `--sample_metadata` | `data/*_wgs_nanopore_wf_metadata.csv` | CSV with cohort-level sample metadata for the `SV_REPORT_INDEX` Quarto report |
| `--knot_config` | `./config_AnnotSV.yaml` | Custom KnotAnnotSV config YAML |
| `--vep_cache` | `data/vep` | Path to pre-built VEP cache (null = skips VEP unless `sv_vcf`/`sv_clinvar` set) |
| `--include_cache` | `false` | Auto-download VEP cache when `vep_cache` is null |
| `--sv_vcf` | `data/gnomad_sv_cnv/gnomad.v4.1.sv.sites.vcf.gz` | gnomAD SV VCF for VEP `--custom` annotation (must be bgzipped + tabix-indexed) |
| `--sv_clinvar` | `data/clinvar/nstd102.GRCh38.variant_call.vcf.gz` | ClinVar SV VCF for VEP `--custom` annotation (must be bgzipped + tabix-indexed) |
| `--sv_reference_filter_policy` | `filter` | What to do when `sv_vcf` / `sv_clinvar` contain non-PASS records: `filter` (annotate against the PASS-only subset), `fail` (abort), `warn` (report only, use as supplied). A per-file census is always written to `reference/sv_filter_check/` |
| `--species` | `homo_sapiens` | VEP species |
| `--cache_version` | `115` | VEP cache version |

> **Note:** `sv_vcf` and `sv_clinvar` must be provided as bgzipped (`.vcf.gz`) files with a co-located tabix index (`.vcf.gz.tbi`). They are checked, and by default reduced to PASS records, before VEP sees them (`--sv_reference_filter_policy`); the raw files on disk are never modified.

---

## Reference data

Reference data comes from two places:

| Kind | Examples | How to get it |
|---|---|---|
| **In the Zenodo archive** | GRCh38 FASTA, gnomAD SV + CNV, ClinVar SV + SNV, EnsDb, CIViC, UCSC chain + cytoband, SAVANA contigs, **and** the project-specific user BEDs, curated gene list and `config_AnnotSV.yaml` | One download — see [Zenodo archive](#zenodo-archive--recommended-starting-point) |
| **Auto-downloadable** | VEP cache, AnnotSV annotations | The pipeline fetches them (see below), or download manually if the network blocks it. Too large for the standard archive |

Between them these cover everything. The archive is the recommended starting point: one download, then the pipeline will fetch VEP and AnnotSV (or stage those manually). Individual download commands are provided if pipeline is unable to retrieve the VEP or AnnotSV references.

Every file, its expected path, and its source URL are listed in [`assets/reference_files_manifest.tsv`](assets/reference_files_manifest.tsv). All paths below are relative to the repository root and match the defaults in `nextflow.config`.

### Automatic downloads

Two references can be fetched by the pipeline itself:

- **AnnotSV annotations** — set `--annotsv_annotations null` and `ANNOTSV_INSTALLANNOTATIONS`
  downloads the bundle (~32 GB) from `www.lbgi.fr` and `ftp.ncbi.nih.gov`.
- **VEP cache** — set `--vep_cache null --include_cache true` and `ENSEMBLVEP_DOWNLOAD`
  fetches the cache (~24 GB) from `ftp.ensembl.org`. With `--include_cache false`,
  VEP instead runs `--offline` using only the `--custom` VCFs.

> **Firewalls and proxies.** Both downloads happen *inside the container*, so the
> host's proxy settings do not apply automatically. On a restricted network set
> `params.proxy_env_file` to an env-file containing `http_proxy`, `https_proxy`,
> `ftp_proxy` and `no_proxy`; it is threaded in as `--env-file` for every process
> that reaches the network (`ENSEMBLVEP_DOWNLOAD`, `ANNOTSV_INSTALLANNOTATIONS`,
> `ANNOTSV_SETUPUSERANNO`, `ANNOTSV_ANNOTSV`, `ENSEMBLVEP_VEP_SNV`, `IGVREPORTS`).


### VEP cache and AnnotSV — manual download

```bash
# ---- VEP cache (~24 GB) -> data/vep/ --------------------------------------
# Release 115 must match `params.cache_version`; bump both together.
mkdir -p data/vep && cd data/vep
wget https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_vep_115_GRCh38.tar.gz          # -> data/vep/homo_sapiens/115_GRCh38/
rm homo_sapiens_vep_115_GRCh38.tar.gz               # optional: the tarball duplicates ~24 GB
cd ../..
# then run with: --vep_cache data/vep

# ---- AnnotSV annotations (~32 GB) -> data/reference/AnnotSV ---------------
# Installs the tool and its human annotation bundle; the pipeline only needs
# the annotations directory ($ANNOTSV/share/AnnotSV).
git clone https://github.com/lgmgeo/AnnotSV.git && cd AnnotSV
make PREFIX=. install && make PREFIX=. install-human-annotation
cd ..
ln -s "$PWD/AnnotSV/share/AnnotSV" data/reference/AnnotSV
# then run with: --annotsv_annotations data/reference/AnnotSV
```

> `data/reference/AnnotSV` may be a **symlink**.

### Additional public downloads

Run from the repository root.

```bash
mkdir -p data/reference data/clinvar data/gnomad_sv_cnv

# ---- Genome FASTA: GRCh38 no-alt analysis set -----------------------------
wget -P data/reference https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# ---- gnomAD SV + CNV (release 4.1) ----------------------------------------
# To update: change "4.1" in the paths. An `aws s3 cp --no-sign-request
# s3://gnomad-public-us-east-1/release/4.1/...` equivalent exists for each file.
wget -P data/gnomad_sv_cnv https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
wget -P data/gnomad_sv_cnv https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
wget -P data/gnomad_sv_cnv https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz

# ---- ClinVar SV (dbVar nstd102) -------------------------------------------
wget -P data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
wget -P data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz.tbi

# ---- ClinVar SNV ----------------------------------------------------------
# Dated release. To update: pick a newer date from the FTP listing and change
# `params.clinvar_snv_vcf` to match — the filename is the version.
wget -P data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260627.vcf.gz
wget -P data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260627.vcf.gz.tbi

# ---- CIViC SNV ------------------------------------------------------------
# CIViC publishes GRCh37 ONLY; PREPARE_REFERENCES lifts it to GRCh38 using the
# UCSC chain below. To update: pick a newer dated release from civicdb.org and
# change `params.civic_snv_vcf`.
wget -P data/reference https://civicdb.org/downloads/01-Jul-2026/01-Jul-2026-civic_accepted.vcf
bgzip data/reference/01-Jul-2026-civic_accepted.vcf
bcftools index -t data/reference/01-Jul-2026-civic_accepted.vcf.gz
wget -P data/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# ---- UCSC cytoband ideogram (IGV reports) ---------------------------------
wget -O data/reference/cytoBandIdeo.hg38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz

# ---- SAVANA contig list ---------------------------------------------------
wget -P data/reference https://raw.githubusercontent.com/cortes-ciriano-lab/savana/refs/heads/main/example/contigs.chr.hg38.txt
```

`EnsDb.Hsapiens.v110.sqlite` (optional, `--ensembl_db`) is not a plain download. It can be generated with Bioconductor AnnotationHub (record `AH113665`) via `ensembldb`.

### Zenodo archive — recommended

| Contents | |
|---|---|
| `data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | 2.9 GB |
| `data/gnomad_sv_cnv/` — gnomAD v4.1 SV + CNV (+ `.tbi`) | 1.6 GB |
| `data/clinvar/` — ClinVar SNV + dbVar nstd102 SV (+ `.tbi`) | 188 MB |
| `data/reference/EnsDb.Hsapiens.v110.sqlite` | 432 MB |
| `data/reference/` — CIViC release, UCSC chain + cytoband, SAVANA contigs | ~400 KB |
| `data/reference/genomic_basis_tall/` — ST16/ST17/ST19 BEDs, ST9 gene list | 1.4 MB |
| `data/reference/HGNC_TCR_genes.tsv`, `config_AnnotSV.yaml` | 64 KB |

The genomic_basis_tall: are **derived from published supplementary tables** (Pölönen et al., *Nature* 632:1082, 2024) and from curated lists, so they cannot be re-downloaded from any URL. 

```bash
# download the archive from the Zenodo record, then, from the repository root:
shasum -a 256 -c SHA256SUMS.curated.*.txt
tar -xzf tall_wgs_refs.curated.<date>.tar.gz
```

To rebuild or refresh it, use [`bin/make_reference_archive.sh`](bin/make_reference_archive.sh). The archive stores **repository-relative paths**, so extracting at the repo root places every file exactly where `nextflow.config` expects it after a user creates the references the first time. 

```bash
bin/make_reference_archive.sh                        # curated tier (default), ~5.2 GB
bin/make_reference_archive.sh --tier curated --dry-run   # list members + sizes only
bin/make_reference_archive.sh --tier bulk            # ~61 GB, adds VEP + AnnotSV
```

`--tier bulk` produces a fully self-contained copy including the VEP cache and
AnnotSV annotations. At ~61 GB it exceeds Zenodo's **50 GB** per-record cap, so it
is split automatically into parts that must be reassembled before extraction:

```bash
cat tall_wgs_refs.bulk.<date>.tar.gz.* > tall_wgs_refs.bulk.<date>.tar.gz
```


---


## Output structure

```
results/
├── merge_sv/
│   ├── bcftools/          # reheader, view, sort, index outputs
│   └── minda/             # MINDA ensemble VCF
├── minda_annotations/
│   ├── minda/             # extracted annotation tables
│   └── bcftools/          # annotated Severus VCF + SAVANA-only VCF
├── annotations_sv/
│   ├── merged_sv_caller/  # final merged+CALLER VCF (Severus + SAVANA-only)
│   ├── annotsv/           # AnnotSV TSV + VCF outputs
│   ├── knotannotsv/       # interactive HTML report
│   └── ensemblvep/        # VEP-annotated VCF (when enabled)
├── annotsv/
│   └── setupuseranno/     # pre-processed AnnotSV reference with user BED annotations
├── filter_sv/             # vembrane-filtered VCF (*.merged.annotated.filtered.pass2.vcf.gz + .tbi)
├── reference/
│   ├── sv_filter_check/   # FILTER census + verdict per SV reference VCF (gnomAD-SV, ClinVar-SV)
│   ├── civic/             # CIViC VCF lifted to GRCh38
│   └── rescue/            # SNV rescue BED
├── cohort/                # per-sample ID-prefixed + region-restricted pre-filter VCFs and record tables (cohort merge inputs)
├── jasmine/               # raw Jasmine merged VCF
├── cohort_merge_sv/       # cohort_sv.regions.bed, cohort.sv.merged.vcf.gz, cohort_sv_recurrent.tsv (+ .map.tsv)
├── genomic_viz/           # inputs/ (links, cna, labels, counts, cohort.maf, cohort.cn.tsv, oncoprint_genes.txt), circos/ (<set>.circos.png|svg), oncoplot/ (cohort.oncoprint.png|svg, onco_matrix.txt, summaries)
├── sv_report/viz/         # the same files re-published beside index.html for the Cohort figures page
├── sv_report/             # cohort-level Quarto HTML report (index.html) + cohort_{sv,cna}_pass2.tsv, cohort_snv.tsv;
│                          #   protein_paint_out/: DISCO HTMLs, sv_tier1_2.svfusion.txt and sv_tier1_2.genefusion.txt (Gene Fusion app input, Format 1)
└── pipeline_info/         # execution timeline, report, trace, DAG
```

---

## Credits

Authored by Jenny L. Smith.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Citations

A list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
