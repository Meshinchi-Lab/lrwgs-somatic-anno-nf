#!/usr/bin/env bash
set -euo pipefail

# Normalize flag-style CSQ in VCF INFO fields to value-style empty-pipe form.
#
# Background: VEP (via the AnnotSV+variantconvert chain) sometimes emits CSQ as
# a value-less flag (";CSQ" instead of ";CSQ=Allele|...|") on breakend records
# whose ALT notation it cannot parse. This is malformed because CSQ is declared
# Type=String, not Flag. vembrane's CSQ parser treats the field as bool=True and
# any expression that touches ANN aborts with:
#
#   Error processing record N: 'bool' object has no attribute 'split'
#
# Earlier "strip" approach: remove the malformed ";CSQ" token entirely. This
# fails silently downstream because vembrane filter (via pysam) re-encodes the
# declared-but-absent CSQ field as flag-style on output, reintroducing the bug.
#
# This "normalize" approach handles TWO failure modes:
#   (a) records that have flag-style ";CSQ" → replace with ";CSQ=|||...|"
#   (b) records that have no CSQ field at all → append ";CSQ=|||...|"
# Both produce records with a syntactically valid value-style CSQ entry that
# matches the subfield count declared in the CSQ header's Format description.
#
# Why both: vembrane's "missing-CSQ" path (triggered by ANN.get() on a record
# whose INFO has no CSQ) writes the record with flag-style ";CSQ" on output.
# So even if STRIP removes all original flag-CSQ, the absent-CSQ records
# generate new flag-CSQ on vembrane filter output — same bug, different cause.
# Pre-populating every record with empty-pipe CSQ closes both holes.
#
# Vembrane and bcftools round-trip the empty-pipe form cleanly, and downstream
# `ANN.get(...)` returns "" for each subfield — which our filter expressions
# handle via `or "0"` / `or ""`.
#
# Implementation: single-pass awk that builds the empty-CSQ replacement when it
# first sees the CSQ header, then applies it to every data record. Single pass
# avoids SIGPIPE under `set -o pipefail` that a header pre-pass would trigger.
#
# Usage: strip_flag_csq.sh <input.vcf[.gz]> <output_prefix>
#   Produces <output_prefix>.vcf.gz (bgzipped).

input="$1"
prefix="$2"

reader=$([[ "$input" == *.gz ]] && echo "zcat" || echo "cat")

$reader "$input" | awk 'BEGIN{FS=OFS="\t"; subfields=0; empty=""}
    /^##INFO=<ID=CSQ,/ {
        pos = index($0, "Format: ")
        if (pos > 0) {
            rest = substr($0, pos + 8)
            end = index(rest, "\"")
            if (end > 1) {
                fmt = substr(rest, 1, end - 1)
                subfields = split(fmt, a, "|")
                empty = "CSQ="
                for (i = 1; i < subfields; i++) empty = empty "|"
            }
        }
        print
        next
    }
    /^#/ { print; next }
    {
        if (subfields == 0) {
            # No CSQ header found; pass through unchanged (no normalization possible)
            print
            next
        }
        # (1) substitute flag-style CSQ → empty-pipe value-style
        sub(/;CSQ$/, ";" empty, $8)
        gsub(/;CSQ;/, ";" empty ";", $8)
        sub(/^CSQ;/, empty ";", $8)
        if ($8 == "CSQ") $8 = empty
        # (2) APPEND empty-pipe CSQ to records that have no CSQ field at all.
        # Vembrane writes a declared-but-absent CSQ as flag-style `;CSQ` when
        # an expression accesses ANN — making every record carry a valid
        # value-style CSQ proactively prevents that re-encoding.
        if ($8 !~ /(^|;)CSQ(=|$)/) {
            if ($8 == "." || $8 == "") {
                $8 = empty
            } else {
                $8 = $8 ";" empty
            }
        }
        print
    }' | bgzip > "${prefix}.vcf.gz"
