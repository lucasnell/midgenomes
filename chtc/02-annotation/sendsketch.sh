#!/bin/bash

#'
#' One-time script (run interactively) using BBmap's sendsketch.sh to
#' look for contamination in the Tgraci and Pstein assemblies.
#'

conda activate main-env


export ASSEMBLY=Tgraci_contigs.fasta.gz
export ONT_READS=basecalls_guppy-5.0.11.fastq.gz

cp /staging/lnell/assemblies/Tgraci_contigs.fasta.gz ./
cp /staging/lnell/ont/basecalls_guppy-5.0.11.fastq.gz ./
cp /staging/lnell/assemblies/Pstein_contigs.fasta.gz ./

mkdir sendsketch

sendsketch.sh in=Tgraci_contigs.fasta.gz out=sendsketch/Tgraci_contigs_sketch.out

sendsketch.sh in=basecalls_guppy-5.0.11.fastq.gz qin=33 reads=100k size=200k \
    out=sendsketch/ont_reads_sketch.out

sendsketch.sh in=Pstein_contigs.fasta.gz out=sendsketch/Pstein_contigs_sketch.out


