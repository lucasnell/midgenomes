#!/bin/bash

#' Use fastp to trim paired-end, RNAseq Illumina reads.
#' See notes below for differences from what you'd do for DNA sequencing.
#'
#' Requires argument for full path to tar file containing reads.
#' Outputs will be `trimmed_${READS_TAR}` and `trimmed_${READS_TAR}.gz`
#' and will be put into the same directory.
#'
#' Usage:
#' fastp-rna.sh READS_TAR
#'


#' The input argument specifies the *.tar file that contains the reads.
#' ** The tar file must be present in `/staging/lnell/ill/rna/` folder.
#' ** It assumes that only reads are inside the tar file and that if the
#'    file names are sorted alphabetically, the first file is #1 of pair,
#'    and the second is #2 of pair.

export READS_TAR_FULL_PATH=$1
if [[ "${READS_TAR_FULL_PATH}" != *.tar ]]; then
    echo -n "ERROR: Argument to fastp-rna.sh must end in *.tar. " 1>&2
    echo "Yours is '${READS_TAR_FULL_PATH}'." 1>&2
    exit 1
fi
if [ ! -f "${READS_TAR_FULL_PATH}" ]; then
    echo "ERROR: '${READS_TAR_FULL_PATH}' does not exist." 1>&2
    exit 1
fi


# For RNAseq, this paper recommends moderate trimming (phred < 5):
# https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full

# Also recommends moderate trimming with read-length filter:
# https://link.springer.com/article/10.1186/s12859-016-0956-2

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate main-env



export READS1=$(read_tar_name ${READS_TAR_FULL_PATH} 1)
export READS2=$(read_tar_name ${READS_TAR_FULL_PATH} 2)

export READS_TAR=$(basename "${READS_TAR_FULL_PATH}")
export READS_TAR_PATH=$(dirname "${READS_TAR_FULL_PATH}")

export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export TRIM_READS_TAR=trimmed_${READS_TAR}
export OUT_DIR=trimmed_${READS_TAR%.tar}

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xf ${READS_TAR_FULL_PATH} -C ./
check_exit_status "extract reads tar file" $?


# The main things happening here are...

# Automatic adapter trimming (on by default)
# polyG tail trimming for NovaSeq sequencing (if relevant, on by default)
# removal of polyA tails (`--trim_poly_x`)
# enable base correction in overlapped regions for PE data (`--correction`)
# no quality filtering (`--disable_quality_filtering`)
# trimming using a low threshold (`--cut_right --cut_right_mean_quality=5`)
# discard too-short reads after trimming (`--length_required=50`)

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --trim_poly_x \
    --correction \
    --disable_quality_filtering \
    --cut_right --cut_right_mean_quality=5 \
    --length_required=50


mkdir ${READS_TAR%.tar}_fastqc

fastqc ${TRIM_READS1} ${TRIM_READS2} -o ${READS_TAR%.tar}_fastqc


tar -cf ${TRIM_READS_TAR} ${TRIM_READS1} ${TRIM_READS2}
mv ${TRIM_READS_TAR} ${READS_TAR_PATH}/


rm ${READS1} ${READS2} ${TRIM_READS1} ${TRIM_READS2}

rm -r '?' 2> /dev/null

cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${READS_TAR_PATH}/

rm -r ${OUT_DIR}



