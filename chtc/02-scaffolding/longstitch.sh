#!/bin/bash


export THREADS=16

. /app/.bashrc

conda activate longstitch-env

# Input FASTA is given by submit file:
export GENOME=$1
if [ ! -f /staging/lnell/${GENOME}.fasta.gz ]; then
    echo -e "\n\nERROR: /staging/lnell/${GENOME}.fasta.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi


# The input file dictates the output name:
export OUT_SUFFIX=L

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME}${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta

# If the output FASTA already exists, this job stops with exit code 0
if [ -f /staging/lnell/${OUT_FASTA}.gz ]; then
    echo -e "\n\nMESSAGE: /staging/lnell/${OUT_FASTA}.gz already exists."
    echo -e "Exiting...\n"
    exit 0
fi


mkdir ${OUT_DIR}

cd ${OUT_DIR}

cp /staging/lnell/${GENOME}.fasta.gz ./ && gunzip ${GENOME}.fasta.gz
# LongStitch requires *.fa ending
mv ${GENOME}.fasta ${GENOME}.fa

export FASTQ=basecalls_guppy-5.0.11
cp /staging/lnell/${FASTQ}.fastq.gz ./
# LongStitch requires *.fq.gz ending
mv ${FASTQ}.fastq.gz ${FASTQ}.fq.gz


# longstitch run
# I chose these values of `k_ntLink` and `w` by trying all combinations of
# `k_ntLink` = 24, 32, 40 and `w` = 100, 250, 500
# These values provided the best combination of good contiguity and
# BUSCO scores.
longstitch tigmint-ntLink-arks \
    draft=${GENOME} reads=${FASTQ} t=${THREADS} G=100000000 \
    k_ntLink=24 w=250

rm -rf ${GENOME}.fa ${FASTQ}.fq.gz

if [ ! -f *-arks.longstitch-scaffolds.fa ]; then
    echo -e "\n\nERROR: LongStitch failed to produce output." 1>&2
    echo -e "Exiting...\n" 1>&2
    cd ..
    rm -r ./${OUT_DIR}
    exit 1
fi


# Moving the final FASTA to staging
# Gets filename of final scaffolds FASTA file
export LS_FASTA=$(basename -- $(readlink -f *-arks.longstitch-scaffolds.fa))
cp $LS_FASTA ${OUT_FASTA}
# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/

rm -r ${OUT_DIR}


