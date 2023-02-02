#!/bin/bash

#' Use fastp to trim paired-end, Poolseq (DNA) Illumina reads.
#'
#' Requires argument for full path to tar file containing reads.
#' Outputs will be `trimmed_${READS_TAR}` and `trimmed_${READS_TAR}.gz`
#' and will be put into the same directory.
#'
#' Usage:
#' fastp-dna.sh READS_TAR
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate main-env


#' The input argument specifies the *.tar file that contains the reads.
#' ** The tar file must be present in `/staging/lnell/ill/dna/raw_fq/` folder.
#' ** It assumes that only reads are inside the tar file and that if the
#'    file names are sorted alphabetically, the first file is #1 of pair,
#'    and the second is #2 of pair.

export READS_TAR=$1
if [[ "${READS_TAR}" != *.tar ]]; then
    echo -n "ERROR: Argument to fastp-dna.sh must end in *.tar. " 1>&2
    echo "Yours is '${READS_TAR}'." 1>&2
    exit 1
fi
if [ ! -f /staging/lnell/ill/dna/raw_fq/${READS_TAR} ]; then
    echo "ERROR: '/staging/lnell/ill/dna/raw_fq/${READS_TAR}' does not exist." 1>&2
    exit 1
fi
export READS1=$(read_tar_name /staging/lnell/ill/dna/raw_fq/${READS_TAR} 1)
export READS2=$(read_tar_name /staging/lnell/ill/dna/raw_fq/${READS_TAR} 2)

export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export TRIM_READS_TAR=trimmed_${READS_TAR}
export OUT_DIR=trimmed_${READS_TAR/.tar/}

mkdir ${OUT_DIR}
cd ${OUT_DIR}


tar -xf /staging/lnell/ill/dna/raw_fq/${READS_TAR} -C ./
check_exit_status "extract reads tar file" $?

# The main things happening here are...
# Automatic adapter trimming (on by default)
# polyG tail trimming for NovaSeq sequencing (if relevant, on by default)
# enable base correction in overlapped regions for PE data (`--correction`)
# no quality filtering (`--disable_quality_filtering`)

# I'm not quality-trimming because `bwa-mem` will soft-mask low-quality reads
# during the alignment phase.

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --thread ${THREADS} \
    --correction \
    --disable_quality_filtering
check_exit_status "fastp" $?

rm ${READS1} ${READS2}

mkdir ${READS_TAR/.tar/}_fastqc

fastqc ${TRIM_READS1} ${TRIM_READS2} -o ${READS_TAR/.tar/}_fastqc

tar -cf ${TRIM_READS_TAR} ${TRIM_READS1} ${TRIM_READS2}
mv ${TRIM_READS_TAR} /staging/lnell/ill/dna/trimmed/

rm ${TRIM_READS1} ${TRIM_READS2}

rm -r '?' 2> /dev/null

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/ill/dna/trimmed/

rm -r ${OUT_DIR}



