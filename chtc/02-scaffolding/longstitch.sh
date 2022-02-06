#!/bin/bash


export THREADS=8

. /app/.bashrc

conda activate longstitch-env

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
    echo "LongStitch failed to produce output." 1>&2
    cd ..
    rm -r ./${OUT_DIR}
    exit 1
fi


# Moving the final FASTA to staging
# Gets filename of final scaffolds FASTA file
export LS_FASTA=$(basename -- $(readlink -f *-arks.longstitch-scaffolds.fa))
cp $LS_FASTA ${OUT_FASTA}
# Keep the uncompressed version for summaries and for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/


# This outputs basics about scaffold sizes:
summ-scaffs.py ${OUT_FASTA}

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS}
conda deactivate


# Now save the whole directory
cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}

