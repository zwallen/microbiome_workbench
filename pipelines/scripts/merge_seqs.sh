#!/bin/bash

# Get inputs
FILE_STRING = $1
OUT_DIR = $2
CPU = $3
MAX_MEM = $4

# Run bbmerge
if echo $FILE_STRING | grep -q ","
then
    ## Extract files
    FILE1 = $(echo $FILE_STRING | awk -F ',' '{print $1}')
    FILE2 = $(echo $FILE_STRING | awk -F ',' '{print $2}')

    ## Derive file name
    FILE_NAME=$(
        echo $FILE1 | \
        awk -F '/' '{print $NF}' | \
        awk -F '_R1' '{print $1}'
    )

    ## Run bbmerge
    bbmerge.sh \
        in1=${FILE1} \
        in2=${FILE2} \
        out=${OUT_DIR}/${FILE_NAME}.fastq.gz \
        adapters=default \
        iterations=5 \
        rem k=62 extend2=50 ecct \
        t=${CPU} -Xmx${MAX_MEM}g \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
else
    ## Give error if two comma separated files not provided
    echo "ERROR: expecting two file paths separated by comma"
    exit 1
fi
