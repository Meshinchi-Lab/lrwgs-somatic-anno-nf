#!/usr/bin/env python3
"""
Merge a vembrane-extracted INFO subset (e.g. CALLER) into the AnnotSV TSV by ID,
drop AnnotSV's bloated full-INFO column, and emit a real `Candidate_genes`
column by intersecting `Gene_name` with a user-supplied gene list.

AnnotSV's `-candidateGenesFile` flag only boosts the `AnnotSV_ranking_score` —
it does NOT emit a per-row list of the candidate genes that hit the SV. We
compute that here so KnotAnnotSV can render a real "Candidate cancer genes"
column instead of the empty placeholder it had before.

AnnotSV emits 1 "full" + N "split" rows per variant, all sharing the same ID.
A left-merge on ID broadcasts the single CALLER value across every row.

Usage: add_caller.py --annotsv-tsv FILE --caller-tsv FILE [--candidate-genes FILE]
                     --output FILE_FULL [--output-html FILE_HTML] [--max-split-rows N]

`--output` writes the FULL merged TSV — fed to KNOTANNOTSV_XL so the Excel
report retains every split-mode gene row. When `--output-html` is given, a
SECOND TSV is written with at most `--max-split-rows` split rows per
AnnotSV_ID (default 10). Split rows whose Gene_name is in the candidate-gene
set are kept first; remaining slots are filled in input order. This file
feeds KNOTANNOTSV_HTML so the interactive browser report stays small enough
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def load_candidate_set(path):
    """Read a one-gene-per-line text file into a set. Strips whitespace and
    skips blank / '#'-prefixed lines. Returns an empty set if path is None,
    empty, or unreadable."""
    if not path:
        return set()
    fp = Path(path)
    if not fp.is_file() or fp.stat().st_size == 0:
        return set()
    return {
        ln.strip()
        for ln in fp.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    }


def filter_splits_per_variant(df, candidate_set, max_split_rows):
    """For each AnnotSV_ID group, keep the full-mode row + at most
    `max_split_rows` split-mode rows. Split rows whose Gene_name is in the
    candidate set are prioritized; remaining slots are filled in input order.
    Variants with ≤ max_split_rows splits pass through unchanged.

    Implementation note: explicit groupby iteration (rather than
    `.apply(...)`) — pandas 2.2+ drops the grouping column from the
    callback's input by default, which strips `AnnotSV_ID` from the output.
    Iterating manually keeps every column intact and is version-stable."""
    if "Annotation_mode" not in df.columns or "AnnotSV_ID" not in df.columns:
        return df

    pieces = []
    for _vid, group in df.groupby("AnnotSV_ID", sort=False):
        splits = group[group["Annotation_mode"] == "split"]
        if len(splits) <= max_split_rows:
            pieces.append(group)
            continue
        non_splits = group[group["Annotation_mode"] != "split"]
        if candidate_set and "Gene_name" in splits.columns:
            cand_mask = splits["Gene_name"].isin(candidate_set)
            ordered = pd.concat([splits[cand_mask], splits[~cand_mask]])
        else:
            ordered = splits
        kept_splits = ordered.head(max_split_rows)
        pieces.append(pd.concat([non_splits, kept_splits]))

    return pd.concat(pieces, ignore_index=True) if pieces else df


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotsv-tsv",     required=True, help="AnnotSV annotated TSV")
    p.add_argument("--caller-tsv",      required=True, help="vembrane TSV with ID + INFO columns to merge")
    p.add_argument("--candidate-genes", default=None,  help="Optional one-gene-per-line file (e.g. OncoKB cancer-gene list). When present, a `Candidate_genes` column is added showing each row's `Gene_name` ∩ this list, and split-mode rows for the HTML output are prioritized by candidate-gene membership.")
    p.add_argument("--output",          required=True, help="Output FULL merged TSV (fed to KNOTANNOTSV_XL).")
    p.add_argument("--output-html",     default=None,  help="Optional second output TSV with at most `--max-split-rows` split rows per AnnotSV_ID (fed to KNOTANNOTSV_HTML).")
    p.add_argument("--max-split-rows",  type=int, default=10, help="Maximum split-mode rows per variant in the HTML output (default 10).")
    args = p.parse_args()

    ann = pd.read_csv(args.annotsv_tsv, sep="\t", dtype=str, keep_default_na=False)
    caller = pd.read_csv(args.caller_tsv, sep="\t", dtype=str, keep_default_na=False)

    if "ID" not in ann.columns:
        sys.exit(f"ERROR: 'ID' column missing from AnnotSV TSV ({args.annotsv_tsv})")
    if "ID" not in caller.columns:
        sys.exit(f"ERROR: 'ID' column missing from caller TSV ({args.caller_tsv})")

    ann = ann.drop(columns=["INFO"], errors="ignore")

    # vembrane's default "long format" prepends a SAMPLE column; AnnotSV already
    # carries Samples_ID, so drop SAMPLE before merging to avoid a redundant column.
    caller = caller.drop(columns=["SAMPLE"], errors="ignore")

    merged = ann.merge(caller, on="ID", how="left")

    # Convert ALT to symbolic notation for non-breakend SVTYPEs. INS records
    # carry the full inserted sequence as ALT (can be many kilobases), which
    # bloats the KnotAnnotSV HTML and is unreadable in the rendered table;
    # DEL/DUP/INV may already be symbolic but are normalized here for
    # consistency. Breakend records — labelled as either `BND` or `TRA`
    # depending on AnnotSV version / caller convention — keep their VCF
    # bracket notation (e.g. `N[chr5:171274715[`) because the partner
    # coordinate is the clinically relevant content and is encoded inline.
    BREAKEND_TYPES = {"BND", "TRA"}
    if "SV_type" in merged.columns and "ALT" in merged.columns:
        non_breakend_mask = ~merged["SV_type"].isin(BREAKEND_TYPES)
        merged.loc[non_breakend_mask, "ALT"] = "<" + merged.loc[non_breakend_mask, "SV_type"] + ">"

    # Compute a real Candidate_genes column. AnnotSV's `-candidateGenesFile`
    # only adjusts the ranking score — it does not emit a per-row list of
    # which candidate genes the SV actually overlaps. Compute that by
    # intersecting the candidate gene set with the *union* of three AnnotSV
    # columns:
    #   - Gene_name    — genes directly overlapping the SV interval
    #   - Closest_left — nearest gene upstream (5' side) if no overlap
    #   - Closest_right— nearest gene downstream (3' side) if no overlap
    # This captures candidate-gene hits for both overlap and proximity —
    # important for BND breakpoints landing in intergenic space near a
    # driver gene (where Gene_name is empty but the fusion partner is
    # still clinically informative).
    #
    # `Gene_name` is populated only on rows where `Annotation_mode == "full"`
    # (split-mode rows have a single gene per row and are aggregated under
    # the same AnnotSV_ID). We compute the intersection on every row,
    # regardless of mode — KnotAnnotSV groups by AnnotSV_ID and shows the
    # full-mode value at the top level.
    candidate_set = load_candidate_set(args.candidate_genes)
    gene_source_cols = [c for c in ("Gene_name", "Closest_left", "Closest_right")
                        if c in merged.columns]
    if gene_source_cols:
        if candidate_set:
            def intersect_genes(row):
                seen = []
                for col in gene_source_cols:
                    cell = row.get(col)
                    if not cell or (isinstance(cell, float) and pd.isna(cell)):
                        continue
                    for g in str(cell).split(";"):
                        g = g.strip()
                        if g and g in candidate_set and g not in seen:
                            seen.append(g)
                return ";".join(seen)
            merged["Candidate_genes"] = merged.apply(intersect_genes, axis=1)
        else:
            # No gene list supplied — emit empty column so the YAML entry
            # still resolves and the KnotAnnotSV layout stays consistent.
            merged["Candidate_genes"] = ""

    # Full output (every split-mode gene row preserved). This feeds the XL
    # report so the spreadsheet retains complete gene-level annotation.
    merged.to_csv(args.output, sep="\t", index=False)

    # HTML-feed output: cap split rows per variant. The interactive browser
    # report is unusable on full-cohort data when a single SV overlaps
    # hundreds of genes (one split row per gene); 60 MB HTML files are
    # routine on real samples. Capping at 10 brings each variant block down
    # to a single screenful while keeping every candidate-gene match in
    # view.
    if args.output_html:
        merged_html = filter_splits_per_variant(merged, candidate_set, args.max_split_rows)
        merged_html.to_csv(args.output_html, sep="\t", index=False)


if __name__ == "__main__":
    main()
