#!/usr/bin/env bash
set -euo pipefail

# Appends the BED file stem to every non-chrom column in each AnnotSV-generated
# .header.tsv, then removes original .bed files from the user annotations tree.
#
# Usage: fix_user_beds.sh <annotations_dir> <genome_build>
#
# AnnotSV's checkBed (AnnotSV-general.tcl) deletes the original .bed file after
# creating the .formatted.sorted.bed, so we iterate over .formatted.sorted.bed
# files rather than the originals. >| overrides noclobber when the header file
# was already present from a prior cp -rL of the annotations directory.

annotations_dir="$1"
genome_build="$2"
user_dir="${annotations_dir}/Annotations_Human/Users/${genome_build}"

for bed_dir in FtIncludedInSV SVincludedInFt AnyOverlap; do
    dir="${user_dir}/${bed_dir}"
    [ -d "$dir" ] || continue
    for fsbed in "$dir"/*.formatted.sorted.bed; do
        [ -f "$fsbed" ] || continue
        stem=$(basename "$fsbed" .formatted.sorted.bed)
        header_file="$dir/${stem}.header.tsv"
        if [ ! -f "$header_file" ]; then
            echo "ERROR: AnnotSV did not create a header file for ${fsbed}" >&2
            echo "ERROR: Check that the original BED file has a '#chrom'-prefixed comment line." >&2
            exit 1
        fi
        comment_line=$(head -1 "$header_file")
        if [ -z "$comment_line" ]; then
            echo "ERROR: Empty header file: ${header_file}" >&2
            echo "ERROR: The original BED file must begin with a tab-delimited comment line." >&2
            exit 1
        fi
        echo "$comment_line" | \
            awk -v stem="$stem" 'BEGIN{FS="\t"; OFS="\t"}{
                out = ""
                for (i=1; i<=NF; i++) {
                    col = $i; sub(/^#/, "", col)
                    if (i == 1) out = "#" col
                    else out = out OFS col "_" stem
                }
                print out
            }' >| "$header_file"
        if [ ! -s "$header_file" ]; then
            echo "ERROR: ${header_file} is empty after stem-appending." >&2
            echo "ERROR: Check that the BED header line uses TAB (not space) delimiters." >&2
            exit 1
        fi
    done
done

# Keep only .formatted.sorted.bed and .header.tsv
find "$user_dir" \
    -name "*.bed" ! -name "*.formatted.sorted.bed" \
    -delete
