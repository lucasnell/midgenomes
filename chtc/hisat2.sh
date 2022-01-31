#!/bin/bash

# Use HISAT2 to align RNAseq reads to assembly to help scaffold

export THREADS=12

. /app/.bashrc
conda activate main-env

# Argument from submit file defines inputs and output names:
export GENOME=$1.fasta
export OUTNAME=rna_hisat2__$1
export BAM=${OUTNAME}.bam


if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


export RNA1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA2=trimmed_TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz


hisat2-build ${GENOME} tany_hisat_idx
hisat2 -x tany_hisat_idx -1 ${RNA1} -2 ${RNA2} -k 3 -p ${THREADS} \
        --pen-noncansplice 1000000 --no-unal | \
    samtools sort -O bam - > ${BAM}

# # Convert to aligned BAM file with index file for storage:
# samtools view -u -h -b ${SAM} | \
#     samtools sort -@ $((THREADS - 2)) -o ${BAM} -
samtools index ${BAM} ${BAM}.bai
mv ${BAM} ${BAM}.bai /staging/lnell/


# (`tany_hisat_idx*` are index files from `hisat2-build`)
rm ${GENOME} ${RNA1} ${RNA2} tany_hisat_idx*


