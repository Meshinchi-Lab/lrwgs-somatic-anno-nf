#!/bin/bash

set -eou pipefail

export NXF_VER=25.10.0
nextflow -v

# INPUT (sample sheet) and OUTDIR are optional overrides; when unset the
# pipeline falls back to `params.input` / `params.outdir` from nextflow.config.
INPUT="${INPUT:-}"
OUTDIR="${OUTDIR:-}"

nextflow \
    -log reports/pipeline.log \
    run main.nf \
    ${INPUT:+--input "$INPUT"} \
    ${OUTDIR:+--outdir "$OUTDIR"} \
    -profile docker,OPBG \
    -cache true \
    -resume \
    -latest
