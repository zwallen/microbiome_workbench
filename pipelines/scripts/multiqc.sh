#!/bin/bash

# Get inputs
FASTQC_DIR = $1
OUT_DIR = $2

# Run multiqc
multiqc \
    $FASTQC_DIR \
    --outdir $OUT_DIR \
    > ${OUT_DIR}/multiqc.log 2>&1