#!/usr/bin/env bash
set -euo pipefail

# Ensures every user BED file in the staged annotations tree has a tab-delimited
# #chrom header line. AnnotSV's checkUsersBED requires this line to name the
# annotation columns; without it, checkUsersBED infers column count from data
# rows and writes an empty .header.tsv, producing unnamed INFO fields in the VCF.
#
# Standard BED6 files (chrom start end name score strand) get the canonical
# header. Files with a different column count get generic col4/col5/... names.
#
# Usage: ensure_bed_headers.sh <annotations_dir> <genome_build>

annotations_dir="$1"
genome_build="$2"
user_dir="${annotations_dir}/Annotations_Human/Users/${genome_build}"

BED6_HEADER=$'#chrom\tchromStart\tchromEnd\tname\tscore\tstrand'

for bed_dir in FtIncludedInSV SVincludedInFt AnyOverlap; do
    dir="${user_dir}/${bed_dir}"
    [ -d "$dir" ] || continue
    for bed in "$dir"/*.bed; do
        [ -f "$bed" ] || continue
        first=$(head -1 "$bed")
        if [[ "$first" == "#chrom"* ]]; then
            continue
        fi
        ncols=$(awk 'NR==1{print NF; exit}' "$bed")
        if [[ "$ncols" -eq 6 ]]; then
            header="$BED6_HEADER"
        else
            header="#chrom"
            for i in $(seq 2 "$ncols"); do
                header="${header}	col${i}"
            done
        fi
        tmp=$(mktemp)
        printf '%s\n' "$header" | cat - "$bed" >| "$tmp"
        mv "$tmp" "$bed"
        echo "Added header to $(basename "$bed") (${ncols} columns)" >&2
    done
done
