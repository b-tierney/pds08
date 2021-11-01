#!/bin/bash

basename=$1
r1=$2
r2=$3

kraken2 --db /athena/masonlab/scratch/users/btt4001/kraken_db --threads 5 --report "$1".kreport --paired $2 $3 > "$1".kraken

bracken -d /athena/masonlab/scratch/users/btt4001/kraken_db -i "$1".kreport-o ${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t ${THRESHOLD}

