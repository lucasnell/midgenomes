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

export READS1=${READ_BASE}_L002_R1_001.fastq.gz
export READS2=${READ_BASE}_L002_R2_001.fastq.gz
export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export OUT_DIR=trimmed_${READ_BASE}


if [ ! -f /staging/lnell/dna/${READS1} ]; then
    echo "/staging/lnell/dna/${READS1} does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/dna/${READS2} ]; then
    echo "/staging/lnell/dna/${READS2} does not exist!" 1>&2
    exit 222
fi


mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/${READS1} ./
cp /staging/lnell/dna/${READS2} ./

# The main things happening here are...
# Automatic adapter trimming (on by default)
# polyG tail trimming for NovaSeq sequencing (on by default)
# enable base correction in overlapped regions for PE data (`--correction`)
# no quality filtering (`--disable_quality_filtering`)

# I'm not quality-trimming because `bwa-mem` will soft-mask low-quality reads
# during the alignment phase.

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --thread ${THREADS} \
    --correction \
    --disable_quality_filtering


mv ${TRIM_READS1} /staging/lnell/dna/trimmed/
mv ${TRIM_READS2} /staging/lnell/dna/trimmed/

rm ${READS1} ${READS2}

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/dna/trimmed/

rm -r ${OUT_DIR}



