#!/usr/bin/env python3
"""Convert bcftools-query TSV of pass2 SVs (and optionally SNVs) into three
bucket-split BEDPE files for igv-reports consumption.

Buckets by inter-anchor distance:
    SMALL:  chrom1 == chrom2 AND end2 - start1 < SMALL_MAX_BP  (default: 1000)
    MEDIUM: chrom1 == chrom2 AND SMALL_MAX_BP <= end2 - start1 < MEDIUM_MAX_BP  (default: 10000)
    LARGE:  chrom1 != chrom2 OR end2 - start1 >= MEDIUM_MAX_BP

SNVs are always routed to SMALL regardless of computed distance.

All string manipulation lives here — no awk or sed in the calling module.

BEDPE Name-column packing
-------------------------
igv-reports v1.16 `BedpeTable` is hardcoded to expose only 6 coordinate
columns + `Name` in the report table; `--info-columns` is not supported for
`.bedpe` input (only `.vcf` / `.bed` / `.maf` / generic tab). To surface
extra annotation in the report, we pack additional pipe-delimited fields
into the `name` column so they appear in the table cell:

    <ID>|<SVTYPE>|<GENES>|caller=<CALLER>|acmg=<ACMG>|rank=<RANK>|st17=<HIT>

Each key=value trailer is omitted when its source value is empty / "." so
compact records don't get a wall of empty fields.

Chromosome-name normalisation
-----------------------------
The SV pass2 VCF and its BND ALT bracket notation use different naming
styles — bare `14` from `%CHROM` versus `chr5` from the ALT bracket parse.
`normalize_chrom()` forces both to the UCSC `chr`-prefix used by BAM,
ST17/ST19 reference BEDs, RefSeq gene tracks, and the cytoband ideogram
so overlap queries and IGV browsing all agree on names.

Candidate-gene-preferred gene list
----------------------------------
When `--candidate-genes` is supplied and one or more of the record's
Gene_name tokens overlap the panel, only the candidate hits render in the
BEDPE Name column. Otherwise the first `GENE_CAP` (10) Gene_name tokens
render with a `; ... (n=<total> genes)` suffix when capped. Matches the
DT-table `summarize_gene_list` convention in `bin/sv_report_index.qmd`.
"""

import argparse
import re
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

# Bucket boundaries for inter-anchor distance (in bp). SVs with same-chrom
# anchors within these ranges route to the corresponding IGVREPORTS_* alias.
SMALL_MAX_BP = 1000    # strict upper bound for SMALL bucket
MEDIUM_MAX_BP = 10000  # strict upper bound for MEDIUM bucket (else -> LARGE)

# Cap on the number of gene symbols rendered in the BEDPE Name column when
# no candidate-gene overlap exists. Matches DT `summarize_gene_list` cap in
# bin/sv_report_index.qmd (user-requested Task 3 convention).
GENE_CAP = 10

# VCF spec BND ALT forms (VCF 4.2 §5.4.1):
#   t[p[  = <ALT>=<REF>[<partner>:<pos>[
#   [p[t  = <ALT>=[<partner>:<pos>[<REF>
#   t]p]  = <ALT>=<REF>]<partner>:<pos>]
#   ]p]t  = <ALT>=]<partner>:<pos>]<REF>
BND_RE = re.compile(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]')


def normalize_chrom(chrom: str) -> str:
    """Force UCSC `chr`-prefix naming on a chromosome label.

    Handles the naming-style split between the SV pass2 VCF (bare `%CHROM`
    like `14`) and its BND ALT bracket notation (already `chr`-prefixed like
    `chr5`). Ensures every BEDPE row has consistent chromosome names that
    match BAM / ST17 / RefSeq / cytoband inputs.

    Rules:
      - already-`chr` prefixed → unchanged
      - `MT` / `mt` / `M` → `chrM` (Ensembl→UCSC mito alias)
      - anything else → prepend `chr`
      - empty / None → returned unchanged
    """
    if not chrom:
        return chrom
    s = str(chrom)
    if s.lower().startswith('chr'):
        return s
    if s.upper() in ('MT', 'M'):
        return 'chrM'
    return 'chr' + s


def parse_bnd_partner(alt: str) -> Tuple[Optional[str], Optional[int]]:
    """Extract (partner_chrom, partner_pos) from a BND ALT string.

    Returns (None, None) if the ALT doesn't contain a bracket-notation
    partner (e.g., symbolic ALT like <DEL> or non-BND ALT)."""
    m = BND_RE.search(alt)
    if not m:
        return (None, None)
    return (m.group(1), int(m.group(2)))


def load_candidate_set(path: Optional[str]) -> Set[str]:
    """Load a one-gene-per-line panel file into a set.

    Blank lines and `#` comments are skipped. Returns an empty set when the
    path is None / empty / missing so downstream code can trivially check
    truthiness. Case-sensitive match against Gene_name tokens (HGNC symbols
    are canonically upper-case)."""
    if not path:
        return set()
    genes: Set[str] = set()
    try:
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                s = line.strip()
                if s and not s.startswith('#'):
                    genes.add(s)
    except OSError:
        # File missing / unreadable — treat as empty panel; script keeps
        # running (Gene_name-first fallback preserves existing behaviour).
        pass
    return genes


def split_gene_tokens(gene: str) -> List[str]:
    """Split an AnnotSV Gene_name cell into a flat unique token list.

    AnnotSV emits full-mode entries `|`-joined and split-mode entries
    comma-appended, e.g. `LOC1|LOC2|LOC3,LOC1,LOC2,LOC3`. Split on both
    delimiters, drop empties / literal `.`, preserve first-seen order,
    de-duplicate."""
    if not gene or gene == '.':
        return []
    parts: List[str] = []
    seen: Set[str] = set()
    # split on either `,` or `|`
    for tok in re.split(r'[,|]', gene):
        tok = tok.strip()
        if tok and tok != '.' and tok not in seen:
            seen.add(tok)
            parts.append(tok)
    return parts


def format_gene_cell(gene: str, candidate_set: Set[str],
                     cap: int = GENE_CAP) -> str:
    """Compose the gene fragment for the BEDPE Name column.

    Precedence (matches DT-table convention in sv_report_index.qmd):
      1. If any Gene_name token overlaps `candidate_set` → return only the
         candidate hits, `; `-joined (no cap — the panel is inherently small).
      2. Else → return first `cap` tokens; when >`cap` tokens exist, append
         ` ... (n=<total> genes)` so the reader sees the truncation depth.
      3. Empty / `.` → `.`.
    """
    tokens = split_gene_tokens(gene)
    if not tokens:
        return '.'
    if candidate_set:
        hits = [t for t in tokens if t in candidate_set]
        if hits:
            return '; '.join(hits)
    if len(tokens) > cap:
        return '; '.join(tokens[:cap]) + f'; ... (n={len(tokens)} genes)'
    return '; '.join(tokens)


def _split_pipe_outside_parens(s: str) -> List[str]:
    """Split on `|` at paren-depth 0.

    The ST17 hit format embeds gene annotations like
    `_(RUNX1|_TLX3::TLX3|_RUNX1)` where the interior `|` is part of one hit,
    not a hit-separator. Mirrors the R `split_pipe_outside_parens()` helper
    in sv_report_index.qmd."""
    if not s:
        return []
    depth = 0
    parts: List[str] = []
    cur = []
    for ch in s:
        if ch == '(':
            depth += 1
            cur.append(ch)
        elif ch == ')':
            depth = max(depth - 1, 0)
            cur.append(ch)
        elif ch == '|' and depth == 0:
            parts.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append(''.join(cur))
    return parts


def _clean_annot_padding(s: str) -> str:
    """Drop AnnotSV `.` / empty tokens produced by full-vs-split padding.

    A cell like `hit1|hit2,hit3,.,.,.,.` becomes `hit1|hit2,hit3`. Mirrors
    the R `clean_annot_padding()` helper. Idempotent — a clean cell round-
    trips unchanged."""
    if not s:
        return ''
    parts = [p.strip() for p in s.split(',')]
    parts = [p for p in parts if p and p != '.']
    return ','.join(parts)


def _pick_top_st17(nm: str, sc: str) -> str:
    """Return the single highest-score `name [score]` pair from a padded
    ST17 name / score cell.

    Steps mirror the qmd `parse_st17_pair()` (top_n=1 flavour for BEDPE
    compactness):
      1. Strip AnnotSV comma padding from both name and score.
      2. Split name via _split_pipe_outside_parens (respects `_(GENE|_...)`).
      3. Split score by top-level `|`.
      4. Pair by index; when scores shorter than names (AnnotSV per-BED-row
         join-level mismatch), repeat the last score.
      5. Return highest-scoring `name [score]`. Ties break by first-seen.
      6. Missing / non-numeric scores fall back to `name` alone.
    Empty string when no meaningful name token exists."""
    nm = _clean_annot_padding(nm)
    sc = _clean_annot_padding(sc)
    if not nm:
        return ''

    names = [n for n in _split_pipe_outside_parens(nm) if n and n != '.']
    if not names:
        return ''

    score_toks = [t.strip() for t in sc.split('|')] if sc else []
    score_toks = [t for t in score_toks if t and t != '.']
    scores: List[Optional[int]] = []
    for t in score_toks:
        try:
            scores.append(int(t))
        except ValueError:
            scores.append(None)
    scores = [s for s in scores if s is not None]

    n_names = len(names)
    if not scores:
        paired: List[Optional[int]] = [None] * n_names
    elif len(scores) >= n_names:
        paired = [scores[i] for i in range(n_names)]
    else:
        # AnnotSV per-BED-row score naturally applies to every gene-pair
        # variant packed inside that BED row's name field.
        paired = list(scores) + [scores[-1]] * (n_names - len(scores))

    # argmax with None-last semantics
    best_i = 0
    best_score = paired[0]
    for i, s in enumerate(paired[1:], start=1):
        if s is None:
            continue
        if best_score is None or s > best_score:
            best_i = i
            best_score = s

    top_name = names[best_i]
    return f'{top_name} [{best_score}]' if best_score is not None else top_name


def pick_st17_hit(svtype: str,
                  name_bnd: str, name_del: str, name_dup: str,
                  name_ins: str, name_inv: str,
                  score_bnd: str, score_del: str, score_dup: str,
                  score_ins: str, score_inv: str) -> str:
    """Pick the SVTYPE-matched ST17 name/score cell and return the top-1
    `name [score]` pair.

    Empty string when no matched hit exists. Handles AnnotSV `,.,.,` padding
    and pipe-inside-parens gene-pair variants via `_pick_top_st17`.
    Mirrors the qmd `pick_st17` + `parse_st17_pair(top_n=1)` construction
    so the BEDPE Name column carries the same recurrence-BED evidence a
    reviewer sees in the DT table."""
    st17_map = {
        'BND': (name_bnd, score_bnd),
        'DEL': (name_del, score_del),
        'DUP': (name_dup, score_dup),
        'INS': (name_ins, score_ins),
        'INV': (name_inv, score_inv),
    }
    nm, sc = st17_map.get(svtype, ('', ''))
    if not nm or nm == '.':
        return ''
    return _pick_top_st17(nm, sc)


def format_name(rid: str, svtype: str, svlen: int, gene_cell: str,
                caller: str, acmg: str, rank: str, st17_hit: str) -> str:
    """Compose the BEDPE Name column with packed annotation fields.

    Format:
        <id>|<svtype>[|+<abs(svlen)>bp for INS]|<genes>[|caller=<X>][|acmg=<Y>][|rank=<Z>][|st17=<W>]

    Optional trailers are skipped when their source value is `.` / empty
    so simple records don't render a wall of blank fields. INS records
    still get the `+Nbp` size marker inserted after SVTYPE (backwards
    compatible with the previous format). See module-docstring rationale."""
    parts: List[str] = [rid, svtype]
    if svtype == 'INS':
        parts.append(f'+{abs(svlen)}bp')
    parts.append(gene_cell)

    def add_kv(key: str, val: str) -> None:
        if val and val not in ('.', 'NA'):
            parts.append(f'{key}={val}')

    add_kv('caller', caller)
    add_kv('acmg',   acmg)
    add_kv('rank',   rank)
    add_kv('st17',   st17_hit)
    return '|'.join(parts)


def classify_bucket(chrom1: str, start1: int, end1: int,
                    chrom2: str, start2: int, end2: int,
                    is_snv: bool) -> str:
    """Route a BEDPE row to 'small', 'medium', or 'large' bucket.

    SNVs (is_snv=True) always land in 'small' regardless of computed distance.
    For non-SNV rows:
        chrom1 != chrom2                          -> 'large'
        chrom1 == chrom2, distance <  1000        -> 'small'
        chrom1 == chrom2, 1000 <= distance < 10000 -> 'medium'
        chrom1 == chrom2, distance >= 10000       -> 'large'
    where distance = max(end1, end2) - min(start1, start2).
    """
    if is_snv:
        return 'small'
    if chrom1 != chrom2:
        return 'large'
    distance = max(end1, end2) - min(start1, start2)
    if distance < SMALL_MAX_BP:
        return 'small'
    if distance < MEDIUM_MAX_BP:
        return 'medium'
    return 'large'


@dataclass
class BedpeRow:
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str
    score: str
    strand1: str
    strand2: str
    svtype: str

    def to_line(self) -> str:
        return '\t'.join(str(x) for x in [
            self.chrom1, self.start1, self.end1,
            self.chrom2, self.start2, self.end2,
            self.name, self.score, self.strand1, self.strand2, self.svtype,
        ]) + '\n'


# Number of tab-separated fields expected in the SV query TSV. See sv_format
# in subworkflows/local/igv_reports/main.nf. First 8 fields (indices 0-7)
# match the pre-refactor layout; fields 8-21 are the new annotation columns.
SV_QUERY_FIELDS = 22


def sv_row_to_bedpe(tsv_row: str, candidate_set: Set[str]) -> Optional[BedpeRow]:
    """Convert one bcftools-query SV row to a BedpeRow.

    Expected column layout (22 fields, tab-separated):
        0  CHROM                (bare, e.g. `14` — normalised to `chr14`)
        1  POS
        2  INFO/END
        3  INFO/SVTYPE
        4  INFO/SVLEN
        5  ALT                  (BND: bracket notation with partner chrom:pos)
        6  ID
        7  INFO/Gene_name       (AnnotSV `|`-joined, split-mode-comma-repeated)
        8  INFO/CALLER
        9  INFO/ACMG_class
       10  INFO/AnnotSV_ranking_score
       11  INFO/CytoBand        (reserved for future use; not packed into Name)
       12  INFO/name_ST17_Alterations.SV.All_BND
       13  INFO/name_ST17_Alterations.SV.All_DEL
       14  INFO/name_ST17_Alterations.SV.All_DUP
       15  INFO/name_ST17_Alterations.SV.All_INS
       16  INFO/name_ST17_Alterations.SV.All_INV
       17  INFO/score_ST17_Alterations.SV.All_BND
       18  INFO/score_ST17_Alterations.SV.All_DEL
       19  INFO/score_ST17_Alterations.SV.All_DUP
       20  INFO/score_ST17_Alterations.SV.All_INS
       21  INFO/score_ST17_Alterations.SV.All_INV

    Returns None (logged) if the BND ALT is unparseable.
    """
    fields = tsv_row.rstrip('\n').split('\t')
    # Older SV VCFs that don't carry the extended annotation set produce
    # short rows — pad with '.' so field access below stays safe.
    if len(fields) < SV_QUERY_FIELDS:
        fields = fields + ['.'] * (SV_QUERY_FIELDS - len(fields))

    (chrom, pos_s, end_s, svtype, svlen_s, alt, rid, gene,
     caller, acmg, rank, _cytoband,
     name_st17_bnd, name_st17_del, name_st17_dup, name_st17_ins, name_st17_inv,
     score_st17_bnd, score_st17_del, score_st17_dup, score_st17_ins, score_st17_inv,
     ) = fields[:SV_QUERY_FIELDS]

    pos = int(pos_s)
    try:
        end = int(end_s) if end_s and end_s != '.' else pos
    except ValueError:
        end = pos
    try:
        svlen_i = int(svlen_s) if svlen_s and svlen_s != '.' else 0
    except ValueError:
        svlen_i = 0

    gene_cell = format_gene_cell(gene, candidate_set)
    st17_hit = pick_st17_hit(
        svtype,
        name_st17_bnd, name_st17_del, name_st17_dup, name_st17_ins, name_st17_inv,
        score_st17_bnd, score_st17_del, score_st17_dup, score_st17_ins, score_st17_inv,
    )
    name = format_name(rid, svtype, svlen_i, gene_cell,
                       caller, acmg, rank, st17_hit)

    chrom_norm = normalize_chrom(chrom)
    if svtype == 'BND':
        p_chrom, p_pos = parse_bnd_partner(alt)
        if p_chrom is None:
            print(f'WARN: unparseable BND ALT for {rid}: {alt!r}', file=sys.stderr)
            return None
        return BedpeRow(
            chrom_norm, pos - 1, pos,
            normalize_chrom(p_chrom), p_pos - 1, p_pos,
            name, '.', '+', '+', 'BND',
        )
    return BedpeRow(
        chrom_norm, pos - 1, pos,
        chrom_norm, end - 1, end,
        name, '.', '+', '+', svtype,
    )


def snv_row_to_bedpe(tsv_row: str, candidate_set: Set[str]) -> BedpeRow:
    """Convert one bcftools-query SNV row (CHROM POS ID CSQ) to a same-anchor
    BedpeRow. VEP CSQ format is pipe-delimited fields with SYMBOL at index 3
    (0-based), and comma-separated across transcripts — we take the first
    transcript's SYMBOL (matches VEP --pick default). CHROM is normalised
    to `chr`-prefix for consistency with SV rows and BAM/BED tracks.

    Gene precedence: if the VEP SYMBOL is in `candidate_set`, keep it as-is
    (single symbol → always shown). SNV records don't have Gene_name multi-
    tokens so the SV cap-at-10 logic doesn't apply."""
    fields = tsv_row.rstrip('\n').split('\t')
    chrom, pos_s, rid, csq = fields[:4]
    pos = int(pos_s)
    gene = '.'
    if csq and csq != '.':
        first_csq = csq.split(',')[0]
        pipe_fields = first_csq.split('|')
        if len(pipe_fields) > 3 and pipe_fields[3]:
            gene = pipe_fields[3]
    name_parts = [rid, 'SNV']
    if gene not in ('.', ''):
        name_parts.append(gene)
        if candidate_set and gene in candidate_set:
            # Explicit marker so a reviewer sees candidate-panel membership
            # at a glance — mirrors SV `st17=` marker convention.
            name_parts.append('candidate=1')
    name = '|'.join(name_parts)
    chrom_norm = normalize_chrom(chrom)
    return BedpeRow(
        chrom_norm, pos - 1, pos,
        chrom_norm, pos - 1, pos,
        name, '.', '+', '+', 'SNV',
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sv-in', required=True,
                    help='bcftools-query TSV: 22 fields (see sv_row_to_bedpe docstring)')
    ap.add_argument('--snv-in',
                    help='Optional bcftools-query TSV for SNVs: CHROM POS ID CSQ')
    ap.add_argument('--candidate-genes',
                    help='Optional one-gene-per-line panel file. When any '
                         'Gene_name token overlaps the panel, only candidate '
                         'hits render in the BEDPE Name column.')
    ap.add_argument('--out-small', required=True)
    ap.add_argument('--out-medium', required=True)
    ap.add_argument('--out-large', required=True)
    args = ap.parse_args()

    candidate_set = load_candidate_set(args.candidate_genes)

    with ExitStack() as stack:
        outs = {
            'small':  stack.enter_context(open(args.out_small,  'w', encoding='utf-8')),
            'medium': stack.enter_context(open(args.out_medium, 'w', encoding='utf-8')),
            'large':  stack.enter_context(open(args.out_large,  'w', encoding='utf-8')),
        }
        with open(args.sv_in, encoding='utf-8') as fin:
            for line in fin:
                if not line.strip():
                    continue
                row = sv_row_to_bedpe(line, candidate_set)
                if row is None:
                    continue
                bucket = classify_bucket(
                    row.chrom1, row.start1, row.end1,
                    row.chrom2, row.start2, row.end2,
                    is_snv=False,
                )
                outs[bucket].write(row.to_line())
        if args.snv_in:
            with open(args.snv_in, encoding='utf-8') as fin:
                for line in fin:
                    if not line.strip():
                        continue
                    row = snv_row_to_bedpe(line, candidate_set)
                    outs['small'].write(row.to_line())
    return 0


if __name__ == '__main__':
    sys.exit(main())
