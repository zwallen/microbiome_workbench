#!/bin/bash

# Get inputs
FILE_STRING = $1
OUT_DIR = $2
CPU = $3
MAX_MEM = $4

# Make directory for filtered sequences
mkdir -p ${OUT_DIR}/BBDuk_Filtered_Sequences

# Run bbduk in paired or single mode depending on input
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

    ## Run bbduk
    bbduk.sh \
        in=${FILE1} \
        in2=${FILE2} \
        out=${OUT_DIR}/${FILE_NAME}_R1.fastq.gz \
        out2=${OUT_DIR}/${FILE_NAME}_R2.fastq.gz \
        outm=${OUT_DIR}/BBDuk_Filtered_Sequences/${FILE_NAME}_R1.fastq.gz \
        outm2=${OUT_DIR}/BBDuk_Filtered_Sequences/${FILE_NAME}_R2.fastq.gz \
        stats=${OUT_DIR}/${FILE_NAME}_stats.txt \
        ftm=5 tpe tbo qtrim=rl trimq=25 minlen=50 ref=adapters,phix \
        entropy=0.5 entropywindow=50 entropyk=5 \
        -Xmx${MAX_MEM}g \
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

    ## Run bbduk
    bbduk.sh \
        in=${FILE} \
        out=${OUT_DIR}/${FILE_NAME}.fastq.gz \
        outm=${OUT_DIR}/BBDuk_Filtered_Sequences/${FILE_NAME}.fastq.gz \
        stats=${OUT_DIR}/${FILE_NAME}_stats.txt \
        ftm=5 qtrim=rl trimq=25 minlen=50 ref=adapters,phix \
        entropy=0.5 entropywindow=50 entropyk=5 \
        -Xmx${MAX_MEM}g \
        > ${OUT_DIR}/${FILE_NAME}.log 2>&1
fi
