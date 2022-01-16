#!/bin/bash

# have job exit if any command returns with non-zero exit status (aka failure)
set -e


export THREADS=8

# Input FASTA is given by submit file:
export GENOME=$1
if [ ! -f /staging/lnell/${GENOME}.fasta.gz ]; then
    echo "/staging/lnell/${GENOME}.fasta.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.fasta.gz ./ && gunzip ${GENOME}.fasta.gz
# LongStitch requires *.fa ending
mv ${GENOME}.fasta ${GENOME}.fa


# The input file dictates the output name:
export OUT_SUFFIX=longstitch

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME}_${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta

mkdir ${OUT_DIR}

mv ${GENOME}.fa ./${OUT_DIR}/

cd ${OUT_DIR}


export FASTQ=basecalls_guppy-5.0.11
cp /staging/lnell/${FASTQ}.fastq.gz ./
# LongStitch requires *.fq.gz ending
mv ${FASTQ}.fastq.gz ${FASTQ}.fq.gz


# longstitch run \
longstitch tigmint-ntLink-arks \
    draft=${GENOME} reads=${FASTQ} t=${THREADS} G=100000000

rm -rf ${GENOME}.fa ${FASTQ}.fq.gz

# Gets filename of final scaffolds FASTA file
export LS_FASTA=$(basename -- $(readlink -f *-arks.longstitch-scaffolds.fa))
cp $LS_FASTA ../

cd ..

# Moving the final FASTA to staging
mv $LS_FASTA ${OUT_FASTA}
./summ_scaffs ${OUT_FASTA}
gzip ${OUT_FASTA}
mv ${OUT_FASTA}.gz /staging/lnell/

# ~~~~~~~~~~~~~
# For now, we'll just delete the main directory. Change this when you
# finalize the pipeline.
# ~~~~~~~~~~~~~
# # Now the whole directory
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR} summ_scaffs

