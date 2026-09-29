#!/usr/bin/env python3
"""
Filter SAVANA VCF to calls absent from MINDA (SAVANA-only calls).

Reads supported_ids.txt produced by minda_extract_annotations.py and removes
all SAVANA records whose ID appears in that set. Also:
  - appends 8 new INFO definitions from header.txt so BCFTOOLS_CONCAT sees
    matching headers in both input VCFs
  - renames the sample column to sample_id for BCFTOOLS_CONCAT compatibility
  - adds INFO/CALLER=SAVANA to each kept record
  - adds INFO/SAVANA_ID=<record.id> for traceability
  - deduplicates same-position INS calls. SAVANA emits multiple INS at the
    same breakpoint with different inserted-sequence lengths; AnnotSV's
    primary key (chrom_start_end_svtype_n) collapses these to one AnnotSV_ID
    and variantconvert then aborts with "Each variant is assumed to only
    have one single line of 'full' annotation" (variantconvert issue #20,
    closed by maintainer as by-design). Per-group selection priority:
    highest QUAL → highest TUMOUR_READ_SUPPORT → longest SVLEN. Missing
    values rank as -inf so any real measurement wins.
"""

import argparse
import gzip
import sys
from pathlib import Path


def open_vcf(path):
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt')
    return open(p, 'r')


def info_get(info_str, key):
    for kv in info_str.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            if k == key:
                return v
        elif kv == key:
            return 'FLAG'
    return None


def selection_key(cols):
    """Sort key for picking the best of several same-position INS records.

    Cascade: highest QUAL → highest TUMOUR_READ_SUPPORT → longest SVLEN.
    Missing / unparseable values rank as -inf so any real measurement wins.
    """
    qual_str = cols[5]
    try:
        qual = float(qual_str) if qual_str not in ('.', '') else float('-inf')
    except ValueError:
        qual = float('-inf')

    info = cols[7]
    try:
        trs = int(info_get(info, 'TUMOUR_READ_SUPPORT') or '-1')
    except ValueError:
        trs = -1
    try:
        svlen = abs(int(info_get(info, 'SVLEN') or '-1'))
    except ValueError:
        svlen = -1

    return (qual, trs, svlen)


def flush_buffer(buf, out):
    """Emit one record per (CHROM, POS): all non-INS pass through, INS deduped.

    Returns (written_count, dropped_ids) for logging.
    """
    if not buf:
        return 0, []
    ins, others = [], []
    for cols in buf:
        if info_get(cols[7], 'SVTYPE') == 'INS':
            ins.append(cols)
        else:
            others.append(cols)
    written = 0
    dropped = []
    for cols in others:
        out.write('\t'.join(cols) + '\n')
        written += 1
    if ins:
        best = max(ins, key=selection_key)
        out.write('\t'.join(best) + '\n')
        written += 1
        for cols in ins:
            if cols is not best:
                dropped.append(cols[2])
    return written, dropped


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--savana_vcf', required=True)
    p.add_argument('--sample_id',  required=True)
    p.add_argument('--supported_ids', help='Required unless --dedup-only is set')
    p.add_argument('--header_txt',    help='Required unless --dedup-only is set')
    p.add_argument('--dedup-only', action='store_true',
                   help='Skip supported_ids filtering and CALLER/SAVANA_ID INFO injection. '
                        'Only dedup same-position INS records. Used by DEDUPLICATE_SAVANA '
                        'subworkflow, which runs before MINDA so supported_ids does not yet exist.')
    args = p.parse_args()

    if not args.dedup_only:
        if not args.supported_ids or not args.header_txt:
            p.error('--supported_ids and --header_txt are required unless --dedup-only is set')
        supported = set(Path(args.supported_ids).read_text().splitlines())
        extra_headers = Path(args.header_txt).read_text().splitlines()
    else:
        supported = set()
        extra_headers = []

    out_suffix = 'savana_dedup' if args.dedup_only else 'savana_only'
    out_path = f"{args.sample_id}.{out_suffix}.vcf"
    written = 0
    all_dropped = []

    # VCF is sorted by (CHROM, POS) so same-position records are consecutive.
    # Buffer records sharing the current key; flush + dedup when the key changes.
    buf = []
    cur_key = None

    with open_vcf(args.savana_vcf) as fh, open(out_path, 'w') as out:
        for line in fh:
            line = line.rstrip('\n')

            if line.startswith('##'):
                out.write(line + '\n')
                continue

            if line.startswith('#CHROM'):
                # Header injection + sample-column rename only apply when feeding
                # BCFTOOLS_CONCAT against the Severus-annotated VCF (filter mode).
                # In dedup-only mode the downstream consumer is MERGE_SV, which
                # expects the original SAVANA header untouched.
                if not args.dedup_only:
                    for h in extra_headers:
                        if h.strip():
                            out.write(h + '\n')
                    cols = line.split('\t')
                    if len(cols) > 9:
                        cols[-1] = args.sample_id
                        line = '\t'.join(cols)
                out.write(line + '\n')
                continue

            # data record
            cols = line.split('\t')
            vid = cols[2]
            if vid in supported:
                continue

            if not args.dedup_only:
                # add CALLER=SAVANA and SAVANA_ID to INFO
                info = cols[7]
                info = (f"CALLER=SAVANA;SAVANA_ID={vid};{info}"
                        if info != '.' else f"CALLER=SAVANA;SAVANA_ID={vid}")
                cols[7] = info

            key = (cols[0], cols[1])
            if key != cur_key:
                w, d = flush_buffer(buf, out)
                written += w
                all_dropped.extend(d)
                buf = []
                cur_key = key
            buf.append(cols)

        # final flush
        w, d = flush_buffer(buf, out)
        written += w
        all_dropped.extend(d)

    label = 'deduped' if args.dedup_only else 'SAVANA-only'
    print(f"Wrote {written} {label} records to {out_path}", file=sys.stderr)
    if all_dropped:
        print(
            f"Deduplicated {len(all_dropped)} same-position INS records "
            f"(dropped IDs): {','.join(all_dropped)}",
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
