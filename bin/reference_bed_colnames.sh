#!/usr/bin/env bash
set -euo pipefail

# Add a standard BED header to every *.bed file in data/reference/genomic_basis_tall/
# that does not already have one. Edits files in-place.

HEADER=$'#chrom\tchromStart\tchromEnd\tname\tscore\tstrand'
BED_DIR="$(dirname "$0")/../data/reference/genomic_basis_tall"

while IFS= read -r -d '' bed; do
    first=$(head -1 "$bed")
    if [[ "$first" == "#chrom"* ]]; then
        echo "skipped (header present): $bed"
        continue
    fi
    tmp=$(mktemp)
    printf '%s\n' "$HEADER" | cat - "$bed" > "$tmp"
    mv "$tmp" "$bed"
    echo "updated: $bed"
done < <(find "$BED_DIR" -maxdepth 1 -name '*.bed' -print0 | sort -z)
