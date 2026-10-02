#!/bin/bash

set -eou pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bin/_nextflow_pinned.sh
. "$SCRIPT_DIR/_nextflow_pinned.sh"

"$NXF_CMD" -v

# INPUT (sample sheet) and OUTDIR are optional overrides; when unset the
# pipeline falls back to `params.input` / `params.outdir` from nextflow.config.
INPUT="${INPUT:-}"
OUTDIR="${OUTDIR:-}"

"$NXF_CMD" \
    -log reports/pipeline.log \
    run main.nf \
    ${INPUT:+--input "$INPUT"} \
    ${OUTDIR:+--outdir "$OUTDIR"} \
    -profile docker,OPBG \
    -cache true \
    -resume \
    -latest
