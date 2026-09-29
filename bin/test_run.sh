#!/bin/bash

set -eou pipefail 

export NXF_VER=25.10.0
nextflow -v 

nextflow \
    -log reports/pipeline.log \
    run main.nf \
    -profile test,docker,emulate_amd64 \
    -with-report reports/pipeline.html \
    -with-dag reports/pipeline.pdf \
    -cache true \
    -resume \
    -latest

