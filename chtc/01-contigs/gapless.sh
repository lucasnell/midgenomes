#!/bin/bash

# Scaffolding and gap-filling using gapless


. /app/.bashrc
conda activate gapless-env
source /staging/lnell/helpers.sh


# Where to send and receive files from/to
export TARGET=/staging/lnell/assemblies

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


# export GENOME=contigs_merged_next_smart.fasta

# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f ${TARGET}/${GENOME}.gz ]; then
    echo -e "\n\nERROR: ${TARGET}/${GENOME}.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi

export OUT_DIR=${GENOME/.fasta/}_gapless
export OUT_FASTA=${OUT_DIR}.fasta


mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp GENOME" $?

export FASTQ=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/${FASTQ} ./
check_exit_status "cp FASTQ" $?




#' gapless.sh [OPTIONS] {long_reads}.fq
#'
#' Parameter	Default	Description
#' -h -?		Display this help and exit
#' -i	(mandatory)	Input assembly (fasta)
#' -j	4	Number of threads
#' -n	3	Number of iterations
#' -o	gapless_run	Output directory (improved assembly is written to gapless.fa in this directory)
#' -r		Restart at the start iteration and overwrite instead of incorporat already present files
#' -s	1	Start iteration (Previous runs must be present in output directory)
#' -t	(mandatory)	Type of long reads (pb_clr,pb_hifi,nanopore)

gapless.sh \
    -i ${GENOME} \
    -j ${THREADS} \
    -o gapless_run \
    -t nanopore \
    ${FASTQ}
check_exit_status "gapless.sh" $?

cp gapless_run/gapless.fa ${OUT_FASTA}


summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads busco

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_DIR} \
    | tee ${OUT_DIR}.csv


# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


