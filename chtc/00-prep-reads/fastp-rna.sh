#!/bin/bash

# Use fastp to trim paired-end, RNAseq Illumina reads.
# See notes below for differences from what you'd do for DNA sequencing

# For RNAseq, this paper recommends moderate trimming (phred < 5):
# https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full

# Also recommends moderate trimming with read-length filter:
# https://link.springer.com/article/10.1186/s12859-016-0956-2

export THREADS=$(count_threads)

. /app/.bashrc
conda activate main-env


#' The input argument specifies the *.tar file that contains the reads.
#' ** The tar file must be present in `/staging/lnell/ill/rna/` folder.
#' ** It assumes that only reads are inside the tar file and that if the
#'    file names are sorted alphabetically, the first file is #1 of pair,
#'    and the second is #2 of pair.

export READS_TAR=$1
if [[ "${READS_TAR}" != *.tar ]]; then
    echo -n "ERROR: Argument to fastp-rna.sh must end in *.tar. " 1>&2
    echo "Yours is '${READS_TAR}'." 1>&2
    exit 1
fi
if [ ! -f /staging/lnell/ill/rna/${READS_TAR} ]; then
    echo "ERROR: '/staging/lnell/ill/rna/${READS_TAR}' does not exist." 1>&2
    exit 1
fi
export READS1=$(read_tar_name /staging/lnell/ill/rna/${READS_TAR} 1)
export READS2=$(read_tar_name /staging/lnell/ill/rna/${READS_TAR} 2)

export TRIM_READS1=trimmed_${READS1}
export TRIM_READS2=trimmed_${READS2}
export TRIM_READS_TAR=trimmed_${READS_TAR}
export OUT_DIR=trimmed_${READS_TAR/.tar/}

mkdir ${OUT_DIR}
cd ${OUT_DIR}

tar -xf /staging/lnell/ill/rna/${READS_TAR} -C ./
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


mkdir ${READS_TAR/.tar/}_fastqc

fastqc ${TRIM_READS1} ${TRIM_READS2} -o ${READS_TAR/.tar/}_fastqc


tar -cf ${TRIM_READS_TAR} ${TRIM_READS1} ${TRIM_READS2}
mv ${TRIM_READS_TAR} /staging/lnell/ill/rna/


rm ${READS1} ${READS2} ${TRIM_READS1} ${TRIM_READS2}

rm -r '?' 2> /dev/null

cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/ill/rna/

rm -r ${OUT_DIR}



