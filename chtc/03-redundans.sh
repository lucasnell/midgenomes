#!/bin/bash

export REDUNDANS=/root/src/redundans/redundans.py
export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
export CONTIGS=midge.contigs.fasta

$REDUNDANS \
        -v -t 12 \
        -l ${LONGREADS} \
        -f ${CONTIGS} \
        -o redundans

