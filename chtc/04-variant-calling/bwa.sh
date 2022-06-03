#!/bin/bash

# export READ_BASE=KS-1-11_S45


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

# Check on BAM file with `bamtools stats`, check status of this call,
# then output the file name
call_bam_stats () {
    local B=$1
    local S=${B/.bam/.stats}
    bamtools stats -in $B | tee $S
    check_exit_status "bamtools stats $2" $?
    echo -e "FILE:" $B "\n**********************************************"
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
export IN_READS_TAR=trimmed_${READ_BASE}.tar
export GENOME=tany_scaffolds.fasta

if [ ! -f /staging/lnell/dna/trimmed/${IN_READS_TAR} ]; then
    echo "/staging/lnell/dna/trimmed/${IN_READS_TAR} does not exist! " 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist! " 1>&2
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
export UNMAPPED_BAM=${READ_BASE}_unmapped.bam
# Intermediates:
export MERGED_READS=merged_trimmed_${READ_BASE}.fastq
export UN_READS1=unmerged_${IN_READS1}
export UN_READS2=unmerged_${IN_READS2}
export MERGED_BAM=${OUT_BAM/.bam/_merged.bam}
export UNMAPPED_MERGED_BAM=${UNMAPPED_BAM/.bam/_merged.bam}
export UN_BAM=${OUT_BAM/.bam/_unmerged.bam}
export UNMAPPED_UN_BAM=${UNMAPPED_BAM/.bam/_unmerged.bam}



#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/trimmed/${IN_READS_TAR} ./ \
    && tar -xf ${IN_READS_TAR} \
    && rm ${IN_READS_TAR}

gunzip ${IN_READS1}.gz && gunzip ${IN_READS2}.gz
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

bwa index ${GENOME}
check_exit_status "bwa-index" $?

# For use later in creating read group info:
export READ_ID=$(head -n 1 ${IN_READS1} | grep -Eo "[ATGCN]+$")

#' ========================================================================
#' Merge paired-end reads.
#' ========================================================================
#'
#' Creating 3 FASTQ files here, for merged reads, unmerged reads # 1,
#' and unmerged reads # 2
bbmerge.sh in1=${IN_READS1} in2=${IN_READS2} out=${MERGED_READS} \
    outu1=${UN_READS1} outu2=${UN_READS2}
check_exit_status "bbmerge" $?

rm ${IN_READS1} ${IN_READS2}


#' ========================================================================
#' Map unmerged reads.
#' ========================================================================
bwa mem -v 1  -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${UN_READS1} ${UN_READS2} | \
    samtools view -@ ${THREADS} -bh - > \
    all_${UN_BAM}
check_exit_status "bwa mem (unmerged)" $?
call_bam_stats all_${UN_BAM} "(unmerged, unfiltered)"

rm ${UN_READS1} ${UN_READS2}

samtools view -@ ${THREADS} -bh -f 0x4 all_${UN_BAM} > ${UNMAPPED_UN_BAM}
check_exit_status "samtools view (unmerged, unmapped)" $?
call_bam_stats ${UNMAPPED_UN_BAM} "(unmerged, unmapped)"

samtools view -@ ${THREADS} -bh -q 20 -F 0x100 all_${UN_BAM} > ${UN_BAM}
check_exit_status "samtools view (unmerged, filtered)" $?
call_bam_stats ${UN_BAM} "(unmerged, filtered)"

rm all_${UN_BAM}


#' ========================================================================
#' Map merged reads.
#' ========================================================================
bwa mem -v 1 -M -t ${THREADS} \
    -R "@RG\tID:${READ_BASE}\tSM:${READ_BASE}_${READ_ID}\tPL:illumina\tLB:lib1" \
    ${GENOME} ${MERGED_READS} | \
    samtools view -@ ${THREADS} -bh - > \
    all_${MERGED_BAM}
check_exit_status "bwa mem (merged)" $?
call_bam_stats all_${MERGED_BAM} "(merged, unfiltered)"

rm ${MERGED_READS}
rm ${GENOME}*

samtools view -@ ${THREADS} -bh -f 0x4 all_${MERGED_BAM} > \
    ${UNMAPPED_MERGED_BAM}
check_exit_status "samtools view (merged, unmapped)" $?
call_bam_stats ${UNMAPPED_MERGED_BAM} "(merged, unmapped)"

samtools view -@ ${THREADS} -bh -q 20 -F 0x100 all_${MERGED_BAM} > ${MERGED_BAM}
check_exit_status "samtools view (merged, filtered)" $?
call_bam_stats ${MERGED_BAM} "(merged, filtered)"

rm all_${MERGED_BAM}





#' ========================================================================
#' Combine BAM files from mappings of merged and unmerged reads.
#' ========================================================================
picard MergeSamFiles \
    I=${MERGED_BAM} \
    I=${UN_BAM} \
    SO=coordinate \
    USE_THREADING=true \
    O=${OUT_BAM} \
    VERBOSITY=WARNING
check_exit_status "Picard_MergeSamFiles" $?

call_bam_stats ${OUT_BAM} "(final)"

rm ${MERGED_BAM} ${UN_BAM}


picard MergeSamFiles \
    I=${UNMAPPED_MERGED_BAM} \
    I=${UNMAPPED_UN_BAM} \
    SO=coordinate \
    USE_THREADING=true \
    O=${UNMAPPED_BAM} \
    VERBOSITY=WARNING
check_exit_status "Picard_MergeSamFiles (unmapped)" $?

call_bam_stats ${UNMAPPED_BAM} "(final, unmapped)"

rm ${UNMAPPED_MERGED_BAM} ${UNMAPPED_UN_BAM}



#' ========================================================================
#' Handle output files
#' ========================================================================

mv ${OUT_BAM} ${TARGET}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
