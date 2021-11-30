#!/bin/bash

# have job exit if any command returns with non-zero exit status (aka failure)
set -e

export OUTDIR=scaffold_longstitch
export THREADS=32


mkdir ${OUTDIR}
cd ${OUTDIR}


export GENOME=haploid_purge_dups
cp /staging/lnell/${GENOME}.fasta.gz ./ && gunzip ${GENOME}.fasta.gz
# LongStitch requires *.fa ending
mv ${GENOME}.fasta ${GENOME}.fa

export FASTQ=basecalls_guppy-5.0.11
cp /staging/lnell/${FASTQ}.fastq.gz ./
# LongStitch requires *.fq.gz ending
mv ${FASTQ}.fastq.gz ${FASTQ}.fq.gz


cp /staging/lnell/ls-env.tar.gz ./

# replace env-name on the right hand side of this line with the name of your
# conda environment
ENVNAME=ls-env
# if you need the environment directory to be named something other than the
# environment name, change this line
ENVDIR=${ENVNAME}

# these lines handle setting up the environment; you shouldn't have to
# modify them
export PATH
mkdir ${ENVDIR}
tar -xzf ${ENVNAME}.tar.gz -C ${ENVDIR}
. ${ENVDIR}/bin/activate


longstitch run draft=${GENOME} reads=${FASTQ} t=${THREADS} G=100000000

rm -rf ${GENOME}.fa ${FASTQ}.fq.gz ${ENVDIR}

cd ..

tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/
rm -r ${OUTDIR}








