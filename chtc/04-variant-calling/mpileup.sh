#!/bin/bash

# export READ_BASE=Lys-19_S14


# Prep mapped reads, then use `samtools mpileup`
# - mark duplicates using `samtools markdup`
# - filter reads (Q >= 20, no secondary alignment) using `samtools view`
# - add group information using `AddOrReplaceReadGroups` in `picard`
# - re-align sequences in the proximity of indels with `IndelRealigner` and
#   `RealignerTargetCreator` in `GATK`
# - `samtools mpileup`

export THREADS=4

# For use as samtools `@` argument
export EXTRA_THREADS=$(($THREADS - 1))

. /app/.bashrc
conda activate main-env


export READ_BASE=$1

export IN_BAM=${READ_BASE}_bwa.bam
export GENOME=tany_scaffolds.fasta

export OUT_DIR=${READ_BASE}_mpileup


if [ ! -f /staging/lnell/dna/bwa/${IN_BAM} ]; then
    echo "/staging/lnell/dna/bwa/${IN_BAM} does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 222
fi


mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/bwa/${IN_BAM} ./
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

# Intermediate output files:

export MARKDUP_OUT=${IN_BAM/.bam/_markdups.bam}
export FILTER_OUT=${MARKDUP_OUT/.bam/_filtered.bam}
export READGROUPS_OUT=${FILTER_OUT/.bam/_readgroups.bam}
export REALIGNED_OUT=${READGROUPS_OUT/.bam/_realigned.bam}


# Indices needed downstream
samtools faidx --length 80 ${GENOME}
picard CreateSequenceDictionary R=${GENOME}

# Mark duplicates
samtools markdup -@ ${EXTRA_THREADS} ${IN_BAM} ${MARKDUP_OUT}

# Filter out reads that are unmapped, secondary alignment, duplicate, or Q < 20
samtools view -h -b -F 4,100,400 -q 20 -@ ${EXTRA_THREADS} \
    -o ${FILTER_OUT} ${MARKDUP_OUT}


picard AddOrReplaceReadGroups \
    I=${FILTER_OUT} \
    O=${READGROUPS_OUT} \
    RGLB=library \
    RGPL=illumina \
    RGPU=barcode \
    RGSM=${READ_BASE}

samtools index -@ ${EXTRA_THREADS} ${READGROUPS_OUT}


GenomeAnalysisTK \
    -T RealignerTargetCreator \
    -nt ${THREADS} \
    -R ${GENOME} \
    -I ${READGROUPS_OUT} \
    -o ${READGROUPS_OUT/.bam/.intervals} \
    -log ${READGROUPS_OUT/.bam/.intervals.log}

# This part can take a while (9+ hours)
GenomeAnalysisTK \
    -T IndelRealigner \
    -R ${GENOME} \
    -I ${READGROUPS_OUT} \
    -targetIntervals ${READGROUPS_OUT/.bam/.intervals} \
    -o ${REALIGNED_OUT}

echo -e "\n\n============================================================\n"
echo -e "Stats on BAM file used in mpileup"
echo -e "\n============================================================\n"
bamtools stats -in ${REALIGNED_OUT} | \
    tee bamtools_stats.out


# Default base quality from
# https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/fq_to_sync_pipeline.sh
# samtools mpileup -B -Q 25 -f ${GENOME} ${REALIGNED_OUT} | \
samtools mpileup -B -f ${GENOME} ${REALIGNED_OUT} | \
    gzip > ${OUT_DIR}.txt.gz

# Just coverage. For an easier time loading into R.
gunzip -c ${OUT_DIR}.txt.gz | cut -f 1,2,4 | \
    gzip > ${OUT_DIR/_mpileup/_coverage}.txt.gz

cp ${OUT_DIR}.txt.gz /staging/lnell/dna/mpileup/
cp ${OUT_DIR/_mpileup/_coverage}.txt.gz /staging/lnell/dna/mpileup/




rm ${IN_BAM} ${GENOME/.fasta/}*

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/dna/mpileup/

rm -r ${OUT_DIR}



