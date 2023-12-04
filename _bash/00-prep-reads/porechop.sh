#!/bin/bash

#' Code to attempt to use PoreChop to remove adapters from P. steinii ONT reads,
#' but it turns out they're already removed.
#'
#' using docker container quay.io/biocontainers/porechop::0.2.4--py39hc16433a_3
#' in an interactive job
#'


export THREADS=$(count_threads)
export IN_READS=Pstein_20180703DL1_LabG03.fastq.gz
export OUT_READS=trimmed_${IN_READS}


cp /staging/lnell/ont/${IN_READS} ./

porechop -i ${IN_READS} \
    --discard_middle --threads ${THREADS} \
    -o ${OUT_READS}

#' Said: "No adapters found - output reads are unchanged from input reads"


