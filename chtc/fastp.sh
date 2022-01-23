#!/bin/bash

# Use fastp to trim paired-end, RNAseq Illumina reads.
# See notes below for differences from what you'd do for DNA sequencing

# For RNAseq, this paper recommends moderate trimming (phred < 5):
# https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full

# Also recommends moderate trimming with read-length filter:
# https://link.springer.com/article/10.1186/s12859-016-0956-2

export THREADS=8


# Argument from submit file gives you the base of the read FASTQ file name.
# Here, I assume that FASTQ file names are of the form
# <sample name>_S<sample number>_L002_R<read number>_001.fastq.gz
# where everything not in brackets is constant.
# This is true for both my RNA and DNA sequencing reads.
# For the file `Ash-19_S5_L002_R2_001.fastq.gz`, the base of the read name
# would be "Ash-19_S5" because everything else can be inferred.

export READ_BASE=$1

export READS1=${READ_BASE}_L002_R1_001.fastq.gz
export READS2=${READ_BASE}_L002_R2_001.fastq.gz
export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export OUT_DIR=trimmed_${READ_BASE}


if [ ! -f /staging/lnell/${READS1} ]; then
    echo "/staging/lnell/${READS1} does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${READS2} ]; then
    echo "/staging/lnell/${READS2} does not exist!" 1>&2
    exit 222
fi


mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/${READS1} ./
cp /staging/lnell/${READS2} ./


# The RNA-specific options here include
# removal of polyA tails (`--trim_poly_x`)
# no quality filtering (`--disable_quality_filtering`)
# trimming using a low threshold (`--cut_right --cut_right_mean_quality=5`)

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --trim_poly_x \
    --length_required=50 \
    --correction \
    --disable_quality_filtering \
    --cut_right --cut_right_mean_quality=5

mv ${TRIM_READS1} /staging/lnell/
mv ${TRIM_READS2} /staging/lnell/

rm ${READS1} ${READS2}


cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/

rm -r ${OUT_DIR}



