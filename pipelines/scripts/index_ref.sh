#!/bin/bash

# Get inputs
REF_FILE = $1
OUT_DIR = $2
CPU = $3
MAX_MEM = $4

# Index reference file with bbsplit
bbsplit.sh \
    ref=${REF_FILE} \
    path=${OUT_DIR} \
    t=${CPU} -Xmx${MAX_MEM}g \
    > ${OUT_DIR}/BBSplit_Indexing.log 2>&1