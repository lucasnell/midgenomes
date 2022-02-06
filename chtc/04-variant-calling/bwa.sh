#!/bin/bash

# Use fastp to trim paired-end, Poolseq (DNA) Illumina reads.


export THREADS=8

. /app/.bashrc
conda activate main-env



# Argument from submit file gives you the base of the read FASTQ file name.
# Here, I assume that FASTQ file names are of the form
# <sample name>_S<sample number>_L002_R<read number>_001.fastq.gz
# where everything not in brackets is constant.
# This is true for both my RNA and DNA sequencing reads.
# For the file `Ash-19_S5_L002_R2_001.fastq.gz`, the base of the read name
# would be "Ash-19_S5" because everything else can be inferred.

export READ_BASE=$1

export OUT_DIR=bwa_${READ_BASE}
export OUT_BAM=${OUT_DIR}.bam

export READS1=trimmed_${READ_BASE}_L002_R1_001.fastq
export READS2=trimmed_${READ_BASE}_L002_R2_001.fastq
export GENOME=

if [ ! -f /staging/lnell/dna/${READS1}.gz ]; then
    echo "/staging/lnell/dna/${READS1}.gz does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/dna/${READS2}.gz ]; then
    echo "/staging/lnell/dna/${READS2}.gz does not exist!" 1>&2
    exit 222
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 222
fi


mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/${READS1}.gz ./ && gunzip ${READS1}.gz
cp /staging/lnell/dna/${READS2}.gz ./ && gunzip ${READS2}.gz
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


# Align, then use `fixmate` for use of samtools markdups later,
# then sort by coordinate:
bwa mem -M -t ${THREADS} ${GENOME} ${READS1} ${READS2} | \
    samtools fixmate -m - - | \
    samtools sort -o sorted.bam > ${OUT_BAM/.bam/_wdups.bam}

# Now mark duplicates
samtools markdup ${OUT_BAM/.bam/_wdups.bam} ${OUT_BAM}

mv ${OUT_BAM} /staging/lnell/dna/


rm ${READS1} ${READS2} ${GENOME} ${OUT_BAM/.bam/_wdups.bam}

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/dna/

rm -r ${OUT_DIR}



