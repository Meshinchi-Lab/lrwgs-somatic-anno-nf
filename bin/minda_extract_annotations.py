#!/usr/bin/env python3
"""
Extract SAVANA annotation data from MINDA ensemble VCF and produce:
  <sample_id>.annotation.tab.gz + .tbi  — tab-delimited, matched by CHROM+POS+~ID
  <sample_id>.header.txt                — 8 ##INFO lines for bcftools annotate
  <sample_id>.savana_supported_ids.txt  — SAVANA IDs seen in MINDA (both BND mates)
"""

import argparse
import gzip
import re
import subprocess
import sys

BND_ALT_RE = re.compile(r'[\[\]]([^:]+):(\d+)[\[\]]')

INFO_LINES = [
    '##INFO=<ID=CALLER,Number=1,Type=String,Description="Originating SV caller">',
    '##INFO=<ID=SAVANA_ID,Number=1,Type=String,Description="SAVANA variant ID from MINDA SUPP_VEC">',
    '##INFO=<ID=SAVANA_SVTYPE,Number=1,Type=String,Description="SAVANA SV type">',
    '##INFO=<ID=SAVANA_POS,Number=1,Type=Integer,Description="SAVANA ALT start position">',
    '##INFO=<ID=SAVANA_END,Number=1,Type=Integer,Description="SAVANA END position">',
    '##INFO=<ID=SAVANA_CHR2,Number=1,Type=String,Description="SAVANA second chromosome (BND only)">',
    '##INFO=<ID=SAVANA_POS2,Number=1,Type=Integer,Description="SAVANA second breakpoint position (BND only)">',
    '##INFO=<ID=SAVANA_MATEID,Number=1,Type=String,Description="SAVANA BND mate ID">',
]


def parse_info(info_str):
    d = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        else:
            d[field] = True
    return d


def open_vcf(path):
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt')
    return open(p, 'r')


def build_savana_dict(savana_vcf):
    """Pass 1: build lookup keyed by SAVANA ID."""
    savana = {}
    with open_vcf(savana_vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos, vid = parts[0], int(parts[1]), parts[2]
            info = parse_info(parts[7])
            svtype = info.get('SVTYPE', '.')
            end = int(info['END']) if 'END' in info else pos
            chr2 = pos2 = mateid = '.'
            if svtype == 'BND':
                m = BND_ALT_RE.search(parts[4])
                if m:
                    chr2, pos2 = m.group(1), int(m.group(2))
                mateid = info.get('MATEID', '.')
            savana[vid] = {
                'chrom': chrom, 'pos': pos, 'svtype': svtype,
                'end': end, 'chr2': chr2, 'pos2': pos2, 'mateid': mateid,
            }
    return savana


def build_severus_dict(severus_vcf):
    """Pass 2: build lookup keyed by SEVERUS ID — ALL SVTYPES."""
    severus = {}
    with open_vcf(severus_vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos, vid = parts[0], int(parts[1]), parts[2]
            info = parse_info(parts[7])
            svtype = info.get('SVTYPE', '.')
            # mate_id only populated for BND; None for DEL/INS/INV/DUP
            mate_id = info.get('MATE_ID', None) if svtype == 'BND' else None
            severus[vid] = {'chrom': chrom, 'pos': pos, 'svtype': svtype, 'mate_id': mate_id}
    return severus


def build_annotation_map(minda_vcf, savana_dict, severus_dict):
    """Pass 3: parse MINDA SUPP_VEC to build annotation_map and supported_ids."""
    annotation_map = {}
    supported_ids = set()

    with open_vcf(minda_vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            info = parse_info(parts[7])
            supp_vec = info.get('SUPP_VEC', '')
            ids = [x for x in supp_vec.split(',') if x]

            severus_id = next((x for x in ids if x.startswith('severus_')), None)
            savana_id  = next((x for x in ids if x.startswith('ID_')), None)

            if not severus_id or not savana_id:
                continue

            sav = savana_dict.get(savana_id)
            if sav is None:
                print(f"WARNING: SAVANA ID {savana_id} not found in SAVANA VCF", file=sys.stderr)
                continue

            ann = {
                'SAVANA_ID':     savana_id,
                'SAVANA_SVTYPE': sav['svtype'],
                'SAVANA_POS':    str(sav['pos']),
                'SAVANA_END':    str(sav['end']),
                'SAVANA_CHR2':   str(sav['chr2']),
                'SAVANA_POS2':   str(sav['pos2']),
                'SAVANA_MATEID': str(sav['mateid']),
            }

            annotation_map[severus_id] = ann
            supported_ids.add(savana_id)

            # propagate mate exclusion for SAVANA BND pairs
            if sav['mateid'] and sav['mateid'] != '.':
                supported_ids.add(sav['mateid'])

            # propagate annotation to SEVERUS BND mate
            sev = severus_dict.get(severus_id)
            if sev and sev['mate_id']:
                annotation_map[sev['mate_id']] = ann

    return annotation_map, supported_ids


def write_annotation_tab(annotation_map, severus_dict, sample_id):
    """Pass 4: write tab-delimited annotation file sorted by chrom/pos, bgzip, tabix.

    bcftools annotate matches this file by CHROM+POS+~ID, bypassing REF/ALT
    allele checks — necessary because Severus uses mixed ALT representations
    (symbolic <DEL>/<INV>, BND bracket notation, and bare dot).
    """
    out_tab = f"{sample_id}.annotation.tab"
    records = []
    for sev_id, sev in severus_dict.items():
        if sev_id in annotation_map:
            ann = annotation_map[sev_id]
            row = (
                sev['chrom'], sev['pos'], sev_id,
                'MINDA',
                ann['SAVANA_ID'], ann['SAVANA_SVTYPE'], ann['SAVANA_POS'],
                ann['SAVANA_END'], ann['SAVANA_CHR2'], ann['SAVANA_POS2'],
                ann['SAVANA_MATEID'],
            )
        else:
            row = (sev['chrom'], sev['pos'], sev_id, 'SEVERUS', '.', '.', '.', '.', '.', '.', '.')
        records.append(row)

    records.sort(key=lambda r: (r[0], r[1]))

    with open(out_tab, 'w') as fh:
        fh.write('#CHROM\tPOS\tID\tCALLER\tSAVANA_ID\tSAVANA_SVTYPE\tSAVANA_POS\t'
                 'SAVANA_END\tSAVANA_CHR2\tSAVANA_POS2\tSAVANA_MATEID\n')
        for row in records:
            fh.write('\t'.join(str(v) for v in row) + '\n')

    subprocess.run(['bgzip', '-f', out_tab], check=True)
    # -s1 CHROM col, -b2 start col, -e2 end col, -c# comment char
    subprocess.run(['tabix', '-s', '1', '-b', '2', '-e', '2', '-c', '#', out_tab + '.gz'], check=True)


def write_header(sample_id):
    """Pass 5a: write header.txt with 8 ##INFO lines."""
    with open(f"{sample_id}.header.txt", 'w') as fh:
        for line in INFO_LINES:
            fh.write(line + '\n')


def write_supported_ids(supported_ids, sample_id):
    """Pass 5b: write savana_supported_ids.txt."""
    with open(f"{sample_id}.savana_supported_ids.txt", 'w') as fh:
        for sid in sorted(supported_ids):
            fh.write(sid + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--minda_vcf',   required=True)
    p.add_argument('--severus_vcf', required=True)
    p.add_argument('--savana_vcf',  required=True)
    p.add_argument('--sample_id',   required=True)
    args = p.parse_args()

    savana_dict  = build_savana_dict(args.savana_vcf)
    severus_dict = build_severus_dict(args.severus_vcf)
    annotation_map, supported_ids = build_annotation_map(
        args.minda_vcf, savana_dict, severus_dict
    )

    write_annotation_tab(annotation_map, severus_dict, args.sample_id)
    write_header(args.sample_id)
    write_supported_ids(supported_ids, args.sample_id)


if __name__ == '__main__':
    main()
