#!/bin/bash

# Get inputs
FILE_STRING = $1
OUT_DIR = $2
CPU = $3
MAX_MEM = $4

# Make directory for filtered sequences
mkdir -p ${OUT_DIR}/BBSplit_Filtered_Sequences

# Run bbsplit in paired or single mode depending on input
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

    ## Run bbsplit
    bbsplit.sh \
        in1=${FILE1} \
        in2=${FILE2} \
        outu1=${OUT_DIR}/${FILE_NAME}_R1.fastq.gz \
        outu2=${OUT_DIR}/${FILE_NAME}_R2.fastq.gz \
        path=${OUT_DIR} \
        basename=${OUT_DIR}/${FILE_NAME}.%_contam_#.fastq.gz \
        refstats=${OUT_DIR}/${FILE_NAME}_refstats.log \
        t=${CPU} -Xmx${MAX_MEM}g \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
else
    ## Copy over file string as file
    FILE = $FILE_STRING

    ## Derive file name
    FILE_NAME=$(
        echo $FILE | \
        awk -F '/' '{print $NF}' | \
        awk -F '.fastq.gz' '{print $1}'
    )

    ## Run bbsplit
    bbsplit.sh \
        in=${FILE} \
        out=${OUT_DIR}/${FILE_NAME}.fastq.gz \
        path=${OUT_DIR} \
        basename=${OUT_DIR}/${FILE_NAME}.%_contam_#.fastq.gz \
        refstats=${OUT_DIR}/${FILE_NAME}_refstats.log \
        t=${CPU} -Xmx${MAX_MEM}g \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
fi

# Move filtered sequences to their own directory
mv ${OUT_DIR}/*contam* ${OUT_DIR}/BBSplit_Filtered_Sequences
