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
#' - Add zero padding to make all contig names the same length
#'
cat ${GENOME/.fasta/_orig.fasta} \
    | tr [:lower:] [:upper:] \
    | tr " " "_" \
    | seqtk rename - contig \
    | seqtk seq -l 60 - \
    | sed "s/^>contig1$/>contig01/g" \
    | sed "s/^>contig2$/>contig02/g" \
    | sed "s/^>contig3$/>contig03/g" \
    | sed "s/^>contig4$/>contig04/g" \
    | sed "s/^>contig5$/>contig05/g" \
    | sed "s/^>contig6$/>contig06/g" \
    | sed "s/^>contig7$/>contig07/g" \
    | sed "s/^>contig8$/>contig08/g" \
    | sed "s/^>contig9$/>contig09/g" \
    > ${GENOME}


gzip ${GENOME}

rm ${GENOME/.fasta/_orig.fasta}
