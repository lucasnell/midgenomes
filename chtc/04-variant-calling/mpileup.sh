#!/bin/bash

#' This script does the following:
#' - mark and remove duplicates using `Picard MarkDuplicates`
#' - re-align sequences in the proximity of indels with `IndelRealigner` and
#'   `RealignerTargetCreator` in `GATK`
#' - `samtools mpileup`


# Check previous command's exit status.
# If != 0, then archive working dir and exit.
check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2" 1>&2
    cd ..
    tar -czf ERROR_${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ERROR_${OUT_DIR}.tar.gz /staging/lnell/dna/mpileup/
    rm -r ${OUT_DIR}
    exit $2
  fi
  echo "Checked step $1"
}

# Print filename (used after `bamtools stats`)
stats_suffix () {
    echo -e "FILE:" $1 "\n**********************************************"
}


export THREADS=4

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


#' ========================================================================
#' Inputs
#' ========================================================================

export IN_BAM=${READ_BASE}_bwa.bam
export GENOME=tany_scaffolds.fasta
if [ ! -f /staging/lnell/dna/bwa/${IN_BAM} ]; then
    echo "/staging/lnell/dna/bwa/${IN_BAM} does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 222
fi


#' ========================================================================
#' Outputs
#' ========================================================================

# Where to send everything when done:
export TARGET=/staging/lnell/dna/mpileup
# Final files / directories
export OUT_DIR=${READ_BASE}_mpileup
export OUT_FILE=${READ_BASE}_mpileup.txt.gz
# Intermediates:
export MARKDUP_OUT=${IN_BAM/.bam/_nodups.bam}
export REALIGNED_OUT=${MARKDUP_OUT/.bam/_realigned.bam}


#' ========================================================================
#' Do the things
#' ========================================================================

#' ----------------------------------------------------
#' Prep for downstream steps.

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/bwa/${IN_BAM} ./
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

# Indices needed downstream
samtools faidx --length 80 ${GENOME}
picard CreateSequenceDictionary R=${GENOME}


#' ----------------------------------------------------
#' Mark and remove duplicates
picard MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=${IN_BAM} \
    O=${MARKDUP_OUT} \
    M=${MARKDUP_OUT/.bam/_report.txt} \
    VALIDATION_STRINGENCY=SILENT
check_exit_status "Picard_MarkDuplicates" $?

rm ${IN_BAM}

bamtools stats -in ${MARKDUP_OUT} | tee ${MARKDUP_OUT/.bam/.stats}
check_exit_status "bamtools stats (nodups)" $?
stats_suffix ${MARKDUP_OUT} | tee ${MARKDUP_OUT/.bam/.stats}

samtools index -@ $(($THREADS - 1)) ${MARKDUP_OUT}
check_exit_status "samtools index (nodups)" $?


#' ----------------------------------------------------
#' Realign around indels

GenomeAnalysisTK \
    -T RealignerTargetCreator \
    -nt ${THREADS} \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -o ${MARKDUP_OUT/.bam/.intervals} \
    -log ${MARKDUP_OUT/.bam/.intervals.log}
check_exit_status "RealignerTargetCreator" $?

# This part can take a while (9+ hours)
GenomeAnalysisTK \
    -T IndelRealigner \
    -R ${GENOME} \
    -I ${MARKDUP_OUT} \
    -targetIntervals ${MARKDUP_OUT/.bam/.intervals} \
    -o ${REALIGNED_OUT}
check_exit_status "IndelRealigner" $?

rm ${MARKDUP_OUT}*

bamtools stats -in ${REALIGNED_OUT} | tee ${REALIGNED_OUT/.bam/.stats}
check_exit_status "bamtools stats (realigned)" $?
stats_suffix ${REALIGNED_OUT} | tee ${REALIGNED_OUT/.bam/.stats}



#' ----------------------------------------------------
#' mpileup
samtools mpileup -B -f ${GENOME} ${REALIGNED_OUT} > ${OUT_FILE/.gz/}
check_exit_status "samtools mpileup" $?

gzip ${OUT_FILE/.gz/}

# Sliding window (500bp with step size of 250bp) of coverage for plotting
window-mpileup.py -r ${GENOME} -s 250 -w 500 ${OUT_FILE}
check_exit_status "window-mpileup.py" $?

# Change default naming to specify window size:
mv ${OUT_FILE/.txt/_window.txt} ${OUT_FILE/.txt/_win500.txt}

rm ${GENOME/.fasta/}*




#' ----------------------------------------------------
#' Handle output files


mv ${OUT_FILE}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

