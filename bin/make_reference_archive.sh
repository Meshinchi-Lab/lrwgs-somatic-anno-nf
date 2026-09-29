#!/usr/bin/env bash
#
# Build Zenodo-ready archives of the pipeline's reference data, preserving the
# repo-relative paths the configs expect, so an archive extracted at the repo
# root drops every file exactly where `nextflow.config` looks for it.
#
# Driven by assets/reference_files_manifest.tsv (basename, config_path, source_url).
#
# TWO TIERS:
#
#   curated  (~5.2 GB)  Everything the pipeline needs EXCEPT the VEP cache and
#                       the AnnotSV annotation bundle. Two kinds of file:
#                         * publicly re-downloadable but awkward/slow -- GRCh38
#                           FASTA, ClinVar dated release, gnomAD SV + CNV,
#                           EnsDb sqlite, CIViC release, UCSC chain + cytoband,
#                           SAVANA contig list;
#                         * project-specific files that exist nowhere else --
#                           the ST9/ST16/ST17/ST19 BEDs derived from the Polonen
#                           2024 supplement, the curated TCR gene list, and the
#                           knotAnnotSV column config.
#                       Fits comfortably in one Zenodo record.
#
#   bulk     (~61 GB)   curated PLUS the VEP cache and AnnotSV annotations, in a
#                       single archive. Both additions are re-downloadable and
#                       already automated by the pipeline (ENSEMBLVEP_DOWNLOAD,
#                       ANNOTSV_INSTALLANNOTATIONS), so this tier exists only for
#                       a fully self-contained copy. It EXCEEDS Zenodo's 50 GB
#                       per-record cap and is therefore split into parts.
#
# Usage:
#   bin/make_reference_archive.sh                        # curated tier (default)
#   bin/make_reference_archive.sh --tier bulk
#   bin/make_reference_archive.sh --tier curated --dry-run
#   bin/make_reference_archive.sh --tier bulk --outdir /scratch --split-size 45G
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$REPO_ROOT/assets/reference_files_manifest.tsv"
OUTDIR="$REPO_ROOT/data/zenodo"   # under data/ which is gitignored
TIER="curated"
SPLIT_SIZE="45G"     # stay under Zenodo's 50 GB per-record cap
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tier)       TIER="$2"; shift 2 ;;
        --outdir)     OUTDIR="$2"; shift 2 ;;
        --split-size) SPLIT_SIZE="$2"; shift 2 ;;
        --manifest)   MANIFEST="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=1; shift ;;
        -h|--help)    sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
case "$TIER" in curated|bulk) ;; *) echo "ERROR: --tier must be curated or bulk" >&2; exit 2 ;; esac

# ── tier membership, keyed on the manifest's config_path column ──────────
# CURATED holds both the awkward-but-public references and the project-specific
# files; BULK_EXTRA is what `--tier bulk` adds on top.
CURATED=(
    # -- project-specific: derived from the Polonen 2024 supplement, not downloadable
    "data/reference/genomic_basis_tall"
    "config_AnnotSV.yaml"
    "data/reference/HGNC_TCR_genes.tsv"
    # -- public, small but version-pinned
    "data/reference/01-Jul-2026-civic_accepted.vcf.gz"
    "data/reference/01-Jul-2026-civic_accepted.vcf.gz.tbi"
    "data/reference/hg19ToHg38.over.chain.gz"
    "data/reference/cytoBandIdeo.hg38.txt.gz"
    "data/reference/contigs.chr.hg38.txt"
    # -- public, large/slow
    "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    "data/clinvar/clinvar_20260627.vcf.gz"
    "data/clinvar/clinvar_20260627.vcf.gz.tbi"
    "data/clinvar/nstd102.GRCh38.variant_call.vcf.gz"
    "data/clinvar/nstd102.GRCh38.variant_call.vcf.gz.tbi"
    "data/gnomad_sv_cnv/gnomad.v4.1.cnv.all.vcf.gz"
    "data/gnomad_sv_cnv/gnomad.v4.1.sv.sites.vcf.gz"
    "data/gnomad_sv_cnv/gnomad.v4.1.sv.sites.vcf.gz.tbi"
    "data/reference/EnsDb.Hsapiens.v110.sqlite"
)
BULK_EXTRA=(
    "data/reference/AnnotSV"      # ~32 GB, symlinked; tar -h dereferences it
    "data/vep/homo_sapiens"       # ~24 GB extracted cache (not the tarball)
)

collect() {
    case "$1" in
        curated) printf '%s\n' "${CURATED[@]}" ;;
        bulk)    printf '%s\n' "${CURATED[@]}" "${BULK_EXTRA[@]}" ;;
    esac
}

mkdir -p "$OUTDIR"
STAMP="$(date +%Y%m%d)"
ARCHIVE="$OUTDIR/tall_wgs_refs.${TIER}.${STAMP}.tar.gz"

# ── resolve members, warn on anything missing rather than failing silently ──
members=(); missing=()
while IFS= read -r rel; do
    [[ -z "$rel" ]] && continue
    if [[ -e "$REPO_ROOT/$rel" ]]; then members+=("$rel"); else missing+=("$rel"); fi
done < <(collect "$TIER")

if (( ${#missing[@]} )); then
    echo "WARNING: ${#missing[@]} path(s) in tier '$TIER' are absent and will be skipped:" >&2
    printf '  - %s\n' "${missing[@]}" >&2
    echo "  Fetch them with the source URL in $(basename "$MANIFEST") before archiving." >&2
fi
(( ${#members[@]} )) || { echo "ERROR: nothing to archive for tier '$TIER'." >&2; exit 1; }

echo "Tier        : $TIER"
echo "Repo root   : $REPO_ROOT"
echo "Archive     : $ARCHIVE"
echo "Split size  : $SPLIT_SIZE (Zenodo per-record cap is 50 GB)"
echo "Members     :"
total_kb=0
for m in "${members[@]}"; do
    # -L: follow symlinks. data/reference/AnnotSV is a symlink OUTSIDE the repo,
    # so an un-dereferenced size (and tar) would be meaningless.
    kb=$(du -skL "$REPO_ROOT/$m" 2>/dev/null | cut -f1); kb=${kb:-0}
    total_kb=$(( total_kb + kb ))
    printf '  %8s  %s\n' "$(du -shL "$REPO_ROOT/$m" 2>/dev/null | cut -f1)" "$m"
done
printf 'Uncompressed total: %s GiB\n' "$(awk -v k="$total_kb" 'BEGIN{printf "%.1f", k/1048576}')"

if (( DRY_RUN )); then echo "[dry-run] stopping before tar."; exit 0; fi

# ── build ────────────────────────────────────────────────────────────────
# -h/--dereference is REQUIRED: data/reference/AnnotSV is a symlink to a path
# outside the repo; without it tar stores a dangling link and the archive is
# useless on any other machine.
# Scratch/OS cruft that lives alongside the curated BEDs but is not reference data.
TAR_EXCLUDES=(
    --exclude='.DS_Store'
    --exclude='temp.bed'
    --exclude='test_user_headers.bed'
    --exclude='*.formatted.sorted.bed'   # regenerated by ANNOTSV_SETUPUSERANNO
    --exclude='.git'
)

echo "Creating archive (this can take a long while for the bulk tier)..."
tar -czhf "$ARCHIVE" "${TAR_EXCLUDES[@]}" -C "$REPO_ROOT" "${members[@]}"

bytes=$(wc -c < "$ARCHIVE")
printf 'Archive size: %s\n' "$(du -h "$ARCHIVE" | cut -f1)"

# ── split if the archive would blow the Zenodo per-record cap ────────────
CAP_BYTES=$(( 50 * 1000 * 1000 * 1000 ))
if (( bytes > CAP_BYTES )); then
    echo "Archive exceeds Zenodo's 50 GB cap — splitting into ${SPLIT_SIZE} parts..."
    # GNU split accepts -d/--additional-suffix; BSD split (macOS) does not and
    # exits with "illegal option". Probe once rather than relying on failure.
    if split --help >/dev/null 2>&1; then
        split -b "$SPLIT_SIZE" -d --additional-suffix=.part "$ARCHIVE" "${ARCHIVE}."
    else
        split -b "$SPLIT_SIZE" "$ARCHIVE" "${ARCHIVE}."        # -> .aa .ab .ac
    fi
    # Both suffix styles sort lexically, so `cat archive.*` reassembles correctly.
    rm -f "$ARCHIVE"
    echo "Parts:"; ls -lh "${ARCHIVE}."* | awk '{print "  "$5"  "$NF}'
    echo "NOTE: parts must be reassembled before extraction:"
    echo "      cat $(basename "$ARCHIVE").* > $(basename "$ARCHIVE")"
fi

# ── checksums + restore instructions ─────────────────────────────────────
( cd "$OUTDIR" && shasum -a 256 "$(basename "$ARCHIVE")"* > "SHA256SUMS.${TIER}.${STAMP}.txt" 2>/dev/null ) || true

cat > "$OUTDIR/RESTORE.md" <<'RESTORE'
# Restoring the reference set

Extract at the **repository root** — the archives store repo-relative paths, so
every file lands where `nextflow.config` expects it.

```bash
cd /path/to/2026-04-30_WGS_Nanopore_Annotation_T-ALL

# verify first
shasum -a 256 -c SHA256SUMS.*.txt

# single-part archive
tar -xzf tall_wgs_refs.curated.<date>.tar.gz

# split archive — reassemble first
cat tall_wgs_refs.bulk.<date>.tar.gz.* > tall_wgs_refs.bulk.<date>.tar.gz
tar -xzf tall_wgs_refs.bulk.<date>.tar.gz
```

## Note on `data/reference/AnnotSV`

In the source checkout this is a **symlink** to a directory outside the repo.
The archive is built with `tar -h` so the real contents travel with it; on
restore it becomes a real directory. That is what the pipeline wants —
`params.annotsv_annotations = "data/reference/AnnotSV"`.

## What is deliberately NOT archived

The VEP cache tarball (`data/vep/homo_sapiens_vep_115_GRCh38.tar.gz`, ~24 GB) is
excluded because the extracted `data/vep/homo_sapiens/` is the same data and is
what VEP actually reads. Re-download it from the URL in
`assets/reference_files_manifest.tsv` only if you want the tarball itself.

Anything in the `mid` or `bulk` tiers can be re-fetched from the `source_url`
column of the manifest instead of being restored from Zenodo.
RESTORE

echo "Wrote $OUTDIR/RESTORE.md and checksums."
