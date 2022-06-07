#!/bin/bash

#'
#' This simplifies the final reference assembly, specifically the contig
#' names and the file name.
#' This can be run on the docker image (run `conda activate main-env`), but
#' it really only needs `seqtk` installed.
#'

export GENOME=tany_contigs.fasta

cp contigs_quickmerge_medaka_nextpolish.fasta.gz ${GENOME/.fasta/_orig.fasta}.gz
gunzip ${GENOME/.fasta/_orig.fasta}.gz

#' Below does the following (in order):
#' - Undo softmasking (I'll be manually softmasking repeat regions later),
#' - Simplify sequence names in headers
#' - Change to narrow (i.e., multi-line) format
#'
cat ${GENOME/.fasta/_orig.fasta} \
    | tr [:lower:] [:upper:] \
    | tr " " "_" \
    | seqtk rename - contig \
    | seqtk seq -l 60 - \
    > ${GENOME}

gzip ${GENOME}

rm ${GENOME/.fasta/_orig.fasta}
