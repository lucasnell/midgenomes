#!/bin/bash

#'
#' Manipulate Illumina reads before using them to polish assembly
#' with NextPolish.
#'

#'
#' This code was run in an interactive job using the following submit file:
#'
# universe =                  docker
# docker_image =              lucasnell/tany_genomics:v0.5.9
# log =                       inter.log
# should_transfer_files =     YES
# when_to_transfer_output =   ON_EXIT
# requirements =              (Target.HasCHTCStaging == true)
# request_cpus =              16
# request_memory =            24GB
# request_disk =              100GB
# queue 1


# Not needed in interactive job:
# . /app/.bashrc

conda activate main-env
source /staging/lnell/helpers.sh


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

mkdir working
cd working


#' Where to send files:
export TARGET=/staging/lnell/dna/trimmed/for-polishing

# Input reads that have already been trimmed for adapters, etc.
export READS1=trimmed_MyKS-19-B_S18_L002_R1_001.fastq.gz
export READS2=trimmed_MyKS-19-B_S18_L002_R2_001.fastq.gz
export READS_TAR=trimmed_MyKS-19-B_S18.tar
cp /staging/lnell/dna/trimmed/${READS_TAR} ./ \
    && tar -xf ${READS_TAR} \
    && rm ${READS_TAR}
check_exit_status "cp, tar, rm reads tar file" $?


#' Filter out reads with any uncalled bases, to prevent N in final assembly.
export NON_READS1=noN_${READS1/.gz/}
export NON_READS2=noN_${READS2/.gz/}
export NON_READS_TAR=noN_${READS_TAR}.gz

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${NON_READS1} --out2 ${NON_READS2} \
    --thread ${THREADS} \
    --disable_length_filtering \
    --disable_adapter_trimming \
    --disable_trim_poly_g \
    --dont_eval_duplication \
    --qualified_quality_phred 0 \
    --unqualified_percent_limit 100 \
    --n_base_limit 0
rm ${READS1} ${READS2} fastp*
tar -czf ${NON_READS_TAR} ${NON_READS1} ${NON_READS2}
mv ${NON_READS_TAR} ${TARGET}/




#' Create another set of reads that are also "normalized" to deal with
#' irregular coverage.
export NORM_READS1=norm_${NON_READS1}
export NORM_READS2=norm_${NON_READS2}
export NORM_READS_TAR=norm_${NON_READS_TAR}

bbnorm.sh in=${NON_READS1} in2=${NON_READS2} \
    out=${NORM_READS1} out2=${NORM_READS2} \
    target=100 mindepth=10 \
    passes=3 \
    -Xmx20g \
    threads=$(( THREADS - 2))
rm ${NON_READS1} ${NON_READS2}
tar -czf ${NORM_READS_TAR} ${NORM_READS1} ${NORM_READS2}
mv ${NORM_READS_TAR} ${TARGET}/


cd ..
rm -r working

