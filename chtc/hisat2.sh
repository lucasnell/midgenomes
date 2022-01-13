#!/bin/bash

# Use HISAT2 to align RNAseq reads to assembly to help scaffold

export THREADS=32


# Argument from submit file defines inputs and output names:
export GENOME=$1.fasta
export OUTNAME=rna_hisat2__$1
export SAM=${OUTNAME}.sam
export BAM=${OUTNAME}.bam


if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


export RNA1=TanyAdult_S1_L002_R1_001.fastq
export RNA2=TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz

# HISAT2 (version 2.2.1)
tar -xzf hisat2-2.2.1.tar.gz
rm hisat2-2.2.1.tar.gz
export PATH=$PATH:$(pwd)/hisat2-2.2.1
# samtools (version 1.14)
tar -xzf samtools-1.14.tar.gz
rm samtools-1.14.tar.gz
export PATH=$PATH:$(pwd)/samtools-1.14/bin


hisat2-build ${GENOME} tany_hisat_idx
hisat2 -x tany_hisat_idx -1 ${RNA1} -2 ${RNA2} -k 3 -p ${THREADS} \
    --pen-noncansplice 1000000 -S ${SAM}

# Convert to aligned BAM file with index file for storage:
samtools view -u -b ${SAM} | \
    samtools sort -@ $((THREADS - 2)) -o ${BAM} -
samtools index ${BAM} ${BAM}.bai
mv ${BAM} ${BAM}.bai /staging/lnell/


# (`tany_hisat_idx*` are index files from `hisat2-build`)
rm -r hisat2-2.2.1 ${RNA1} ${RNA2} ${GENOME} ${SAM} tany_hisat_idx*

