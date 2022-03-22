#!/bin/bash


#' Align Pool-seq reads to genome assembly using `bwa mem`.
#'
#' It first merges paired-end reads, then separately aligns both merged reads
#' and any reads that can't be merged.
#' It then combines these alignments into a single sorted BAM file
#' that also contains read group info.

# Check previous command's exit status.
# If != 0, then archive working dir and exit.
check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2" 1>&2
    cd ..
    tar -czf ERROR_${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ERROR_${OUT_DIR}.tar.gz /staging/lnell/dna/bwa/
    rm -r ${OUT_DIR}
    exit $2
  fi
  echo "Checked step $1"
}

# Print filename (used after `bamtools stats`)
stats_suffix () {
    echo -e "FILE:" $1 "\n**********************************************"
}


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


#' ========================================================================
#' Inputs
#' ========================================================================

export IN_READS1=trimmed_${READ_BASE}_L002_R1_001.fastq
export IN_READS2=trimmed_${READ_BASE}_L002_R2_001.fastq
export GENOME=tany_scaffolds.fasta

if [ ! -f /staging/lnell/dna/trimmed/${IN_READS1}.gz ]; then
    echo "/staging/lnell/dna/trimmed/${IN_READS1}.gz does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/dna/trimmed/${IN_READS2}.gz ]; then
    echo "/staging/lnell/dna/trimmed/${IN_READS2}.gz does not exist!" 1>&2
    exit 222
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 222
fi


#' ========================================================================
#' Outputs
#' ========================================================================

# Where to send everything when done:
export TARGET=/staging/lnell/dna/bwa
# Final files / directories
export OUT_DIR=${READ_BASE}_bwa
export OUT_BAM=${READ_BASE}_bwa.bam
# Intermediates:
export MERGED_READS=merged_trimmed_${READ_BASE}.fastq
export UN_READS1=unmerged_${IN_READS1}
export UN_READS2=unmerged_${IN_READS2}
export MERGED_BAM=${OUT_BAM/.bam/_merged.bam}
export UN_BAM=${OUT_BAM/.bam/_unmerged.bam}


#' ========================================================================
#' Do the things
#' ========================================================================

#' ----------------------------------------------------
#' Prep for downstream steps.

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/trimmed/${IN_READS1}.gz ./ && gunzip ${IN_READS1}.gz
cp /staging/lnell/dna/trimmed/${IN_READS2}.gz ./ && gunzip ${IN_READS2}.gz
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

bwa index ${GENOME}
check_exit_status "bwa-index" $?

# For use later in creating read group info:
export READ_ID=$(head -n 1 ${IN_READS1} | grep -Eo "[ATGCN]+$")

#' ----------------------------------------------------
#' Merge paired-end reads.
#' Creating 3 FASTQ files here, for merged reads, unmerged reads # 1,
#' and unmerged reads # 2
bbmerge.sh in1=${IN_READS1} in2=${IN_READS2} out=${MERGED_READS} \
    outu1=${UN_READS1} outu2=${UN_READS2}
check_exit_status "bbmerge" $?

rm ${IN_READS1} ${IN_READS2}


#' ----------------------------------------------------
#' Map unmerged reads.
bwa mem -v 1  -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${UN_READS1} ${UN_READS2} | \
    samtools view -@ ${THREADS} -bh -q 20 -F 0x100 - > \
    ${UN_BAM}
check_exit_status "bwa mem (unmerged)" $?

bamtools stats -in ${UN_BAM} | tee ${UN_BAM/.bam/.stats}
check_exit_status "bamtools stats (unmerged)" $?
stats_suffix ${UN_BAM} | tee ${UN_BAM/.bam/.stats}

rm ${UN_READS1} ${UN_READS2}


#' ----------------------------------------------------
#' Map merged reads.
bwa mem -v 1  -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${MERGED_READS} | \
    samtools view -@ ${THREADS} -bh -q 20 -F 0x100 - > \
    ${MERGED_BAM}
check_exit_status "bwa mem (merged)" $?

bamtools stats -in ${MERGED_BAM} | tee ${MERGED_BAM/.bam/.stats}
check_exit_status "bamtools stats (merged)" $?
stats_suffix ${MERGED_BAM} | tee ${MERGED_BAM/.bam/.stats}

rm ${MERGED_READS}
rm ${GENOME}*


#' ----------------------------------------------------
#' Combine BAM files from mappings of merged and unmerged reads.
picard MergeSamFiles \
    I=${MERGED_BAM} \
    I=${UN_BAM} \
    SO=coordinate \
    USE_THREADING=true \
    O=${OUT_BAM}
check_exit_status "Picard_MergeSamFiles" $?

bamtools stats -in ${OUT_BAM} | tee ${OUT_BAM/.bam/.stats}
check_exit_status "bamtools stats (final)" $?
stats_suffix ${OUT_BAM} | tee ${OUT_BAM/.bam/.stats}

rm ${MERGED_BAM} ${UN_BAM}


#' ----------------------------------------------------
#' Handle output files

mv ${OUT_BAM} ${TARGET}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
