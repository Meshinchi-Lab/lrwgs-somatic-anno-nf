#!/bin/bash

set -eou pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bin/_nextflow_pinned.sh
. "$SCRIPT_DIR/_nextflow_pinned.sh"

"$NXF_CMD" -v

"$NXF_CMD" \
    -log reports/pipeline.log \
    run main.nf \
    -profile test,docker,emulate_amd64 \
    -with-report reports/pipeline.html \
    -with-dag reports/pipeline.pdf \
    -cache true \
    -resume \
    -latest
