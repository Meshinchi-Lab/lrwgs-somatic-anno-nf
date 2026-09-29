#!/usr/bin/env bash
#
# Download Oxford Nanopore WGS reads for the REH B-ALL cell line from ENA.
#
#   Study : PRJNA600820  "REH MultiOmics"  (Uppsala University)
#   Sample: SAMN13831871
#   Use   : open access, no data-use agreement, no registration.
#
# WHY NO POD5/FAST5 --------------------------------------------------------
# Raw signal is NOT available for this study, so re-basecalling is impossible.
# Verified 2026-09-16 against both archives:
#   * ENA  `submitted_ftp` / `submitted_format` are EMPTY for all 9 ONT runs.
#     Of 33 runs in the study only 4 carry submitted files, and all 4 are
#     ILLUMINA BAM;BAI (10x phased_possorted / GemCode / ATAC / recal BAM).
#   * NCBI SRA lists exactly two files for SRR22730978 -- "SRA Normalized" and
#     the fastq `REH_WGS_ONT_PromethION_10kb.pass.fastq.gz`. The string "pod5"
#     and "fast5" appear nowhere in the SRA XML.
# The submitted filename also shows the reads are already "pass"-filtered from
# a ~10 kb library, i.e. basecalled and quality-filtered by the submitter.
#
# Re-check before assuming this is still true:
#   curl -s -G "https://www.ebi.ac.uk/ena/portal/api/search" \
#     --data-urlencode 'result=read_run' \
#     --data-urlencode 'query=run_accession="SRR22730978"' \
#     --data-urlencode 'fields=submitted_format,submitted_ftp' \
#     --data-urlencode 'format=tsv'
#
# Usage:
#   bin/download_reh_ont.sh                      # deepest run only (~84 GB)
#   bin/download_reh_ont.sh --run SRR22730978    # a specific run
#   bin/download_reh_ont.sh --all                # all 9 ONT WGS runs (~271 GB)
#   bin/download_reh_ont.sh --list               # show runs, sizes, exit
#   bin/download_reh_ont.sh --outdir /scratch/reh --verify-only
#
set -euo pipefail

OUTDIR="data/reh_ont"
RUNS=""
DO_ALL=0
LIST_ONLY=0
VERIFY_ONLY=0

# run_accession | ~coverage | fastq bytes | md5 | ftp path
MANIFEST="\
SRR22730978|31x|84148319205|91f21daea6d199bcfcbe4ce79999cc4d|ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/078/SRR22730978/SRR22730978_1.fastq.gz
SRR21147769|18x|54770367775|7e48c42ab1f3702904511bef0da0c45b|ftp.sra.ebi.ac.uk/vol1/fastq/SRR211/069/SRR21147769/SRR21147769_1.fastq.gz
SRR23704822|18x|54428344803|e8fdf2f5e5f6a4b3cb8489e66beae643|ftp.sra.ebi.ac.uk/vol1/fastq/SRR237/022/SRR23704822/SRR23704822_1.fastq.gz
SRR23054498|18x|48963843296|b902013819796e27aae12c172b2e6e11|ftp.sra.ebi.ac.uk/vol1/fastq/SRR230/098/SRR23054498/SRR23054498_1.fastq.gz
SRR22444743|5x|16876027837|35f275bbc900333ab8e48e0fb7791aea|ftp.sra.ebi.ac.uk/vol1/fastq/SRR224/043/SRR22444743/SRR22444743_1.fastq.gz
SRR22444742|2x|6457014778|79a0683396596315a47a7014cee24a63|ftp.sra.ebi.ac.uk/vol1/fastq/SRR224/042/SRR22444742/SRR22444742_1.fastq.gz
SRR23704826|1x|887230664|5efd21b736f88fa35deab16e810de2d7|ftp.sra.ebi.ac.uk/vol1/fastq/SRR237/026/SRR23704826/SRR23704826_1.fastq.gz
SRR23704825|1x|867414985|3203b1b1fa364ce96cb67846f893fb4e|ftp.sra.ebi.ac.uk/vol1/fastq/SRR237/025/SRR23704825/SRR23704825_1.fastq.gz
SRR22444744|1x|3205090323|0e92fa314d074efc488d156e2d22e6d1|ftp.sra.ebi.ac.uk/vol1/fastq/SRR224/044/SRR22444744/SRR22444744_1.fastq.gz"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)         [[ $# -ge 2 ]] || { echo "ERROR: --run needs a value" >&2; exit 2; }
                       RUNS="$RUNS $2"; shift 2 ;;
        --outdir)      [[ $# -ge 2 ]] || { echo "ERROR: --outdir needs a value" >&2; exit 2; }
                       OUTDIR="$2"; shift 2 ;;
        --all)         DO_ALL=1; shift ;;
        --list)        LIST_ONLY=1; shift ;;
        --verify-only) VERIFY_ONLY=1; shift ;;
        -h|--help)     sed -n '2,36p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if (( LIST_ONLY )); then
    printf "  %-13s %-5s %10s  %s\n" RUN COV SIZE FILE
    while IFS='|' read -r acc cov bytes md5 url; do
        printf "  %-13s %-5s %7.1f GB  %s\n" "$acc" "$cov" "$(echo "$bytes/1000000000" | bc -l)" "${url##*/}"
    done <<< "$MANIFEST"
    echo
    echo "  POD5/FAST5: not deposited for this study (see header)."
    exit 0
fi

(( DO_ALL )) && RUNS=$(awk -F'|' '{print $1}' <<< "$MANIFEST")
[[ -z "${RUNS// }" ]] && RUNS="SRR22730978"      # default: deepest run

mkdir -p "$OUTDIR"

# Prefer curl; fall back to wget. Both resume partial transfers, which matters
# a great deal for an 84 GB file over FTP.
if command -v curl >/dev/null 2>&1; then DL=curl
elif command -v wget >/dev/null 2>&1; then DL=wget
else echo "ERROR: need curl or wget" >&2; exit 1; fi

# GNU coreutils vs BSD: md5sum on Linux, md5 -q on macOS.
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    else md5 -q "$1"; fi
}

total_fail=0
for acc in $RUNS; do
    line=$(grep "^${acc}|" <<< "$MANIFEST" || true)
    if [[ -z "$line" ]]; then
        echo "WARNING: $acc is not an ONT WGS run in PRJNA600820 — skipping" >&2
        continue
    fi
    IFS='|' read -r _ cov bytes md5 url <<< "$line"
    dest="$OUTDIR/${url##*/}"

    printf '\n=== %s (~%s, %.1f GB) ===\n' "$acc" "$cov" "$(echo "$bytes/1000000000" | bc -l)"

    if (( ! VERIFY_ONLY )); then
        if [[ -f "$dest" ]] && [[ "$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest")" == "$bytes" ]]; then
            echo "  already complete — skipping download"
        else
            echo "  downloading from ENA (resumable)..."
            case "$DL" in
                # -C - resumes; --retry survives transient FTP drops
                curl) curl -L -C - --retry 5 --retry-delay 10 --fail -o "$dest" "https://$url" ;;
                wget) wget -c --tries=5 --waitretry=10 -O "$dest" "https://$url" ;;
            esac
        fi
    fi

    [[ -f "$dest" ]] || { echo "  ERROR: $dest not present" >&2; total_fail=1; continue; }

    echo "  verifying md5 (this reads the whole file)..."
    got=$(md5_of "$dest")
    if [[ "$got" == "$md5" ]]; then
        echo "  md5 OK  $got"
    else
        echo "  md5 MISMATCH: expected $md5, got $got" >&2
        echo "  delete the file and re-run to re-download." >&2
        total_fail=1
    fi
done

echo
if (( total_fail )); then
    echo "Finished WITH ERRORS — see messages above." >&2
    exit 1
fi
cat <<'NEXT'
Done. Next steps (reads are already basecalled; no POD5 to re-call):

  minimap2 -ax map-ont -t 16 --MD -Y \
      data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      <fastq.gz> | samtools sort -@ 8 -o REH.ont.bam -
  samtools index REH.ont.bam

REH carries the ETV6-RUNX1 t(12;21) fusion — a useful positive control for the
SV flow. The same study also has PacBio WGS, Illumina WGS and RNA-seq on this
cell line, which can serve as an orthogonal truth set.
NEXT
