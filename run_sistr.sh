#!/bin/bash

PREVIEW="${1:-}"

nextflow run sistr_pipeline.nf \
    --input_s3 "s3://your-bucket-name/validationset" \
    --outdir "s3://your-bucket-name/sistr_results" \
    $PREVIEW