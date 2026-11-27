#!/bin/bash

# Get inputs
FILE_STRING=$1
OUT_DIR=$2

# Run fastqc in paired or single mode depending on input
if echo $FILE_STRING | grep -q ","
then
    ## Extract files
    FILE1=$(echo $FILE_STRING | awk -F ',' '{print $1}')
    FILE2=$(echo $FILE_STRING | awk -F ',' '{print $2}')

    ## Derive file name
    FILE_NAME=$(
        echo $FILE1 | \
        awk -F '/' '{print $NF}' | \
        awk -F '_R1' '{print $1}'
    )

    ## Run fastqc
    fastqc \
        $FILE1 $FILE2 \
        -d $OUT_DIR \
        --outdir $OUT_DIR \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
else
    ## Copy over file string as file
    FILE=$FILE_STRING

    ## Derive file name
    FILE_NAME=$(
        echo $FILE | \
        awk -F '/' '{print $NF}' | \
        awk -F '.fastq.gz' '{print $1}'
    )

    ## Run fastqc
    fastqc \
        $FILE \
        -d $OUT_DIR \
        --outdir $OUT_DIR \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
fi
