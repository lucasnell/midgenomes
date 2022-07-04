#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using RNAseq alignments (pipeline B)
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation

export OUT_DIR=tany_braker_rna

mkdir working
cd working

export GENOME=tany_contigs_masked.fasta

cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp genome" $?


#' ===========================================================================
#' ===========================================================================
#'
#' Align RNAseq reads to assembly
#'
#' ===========================================================================
#' ===========================================================================

conda activate main-env


hisat2-build -q ${GENOME} tany_hisat_idx
check_exit_status "hisat2-build" $?


#' ----------------------
#' Adults

export RNA_READS_ADULT1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA_READS_ADULT2=trimmed_TanyAdult_S1_L002_R2_001.fastq
export RNA_READS_ADULT_TAR=trimmed_TanyAdult_S1.tar
export BAM_ADULT=tany_rna_adults.bam

tar -xf /staging/lnell/ill/rna/${RNA_READS_ADULT_TAR} -C ./
check_exit_status "extract rna adult" $?

gunzip ${RNA_READS_ADULT1}.gz \
    && gunzip ${RNA_READS_ADULT2}.gz
check_exit_status "gunzip rna - adults" $?

# (Pipe to samtools is to convert it to an aligned BAM file)
hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_ADULT1} -2 ${RNA_READS_ADULT2} \
    -k 3 -p ${THREADS} \
    | samtools sort -O bam - \
    > ${BAM_ADULT}
check_exit_status "hisat2 - adults" $?

cp ${BAM_ADULT} ${TARGET}/

rm ${RNA_READS_ADULT1} ${RNA_READS_ADULT2}


#' ----------------------
#' Juveniles

export RNA_READS_JUVEN1=trimmed_TanyJuven_S2_L002_R1_001.fastq
export RNA_READS_JUVEN2=trimmed_TanyJuven_S2_L002_R2_001.fastq
export RNA_READS_JUVEN_TAR=trimmed_TanyJuven_S2.tar
export BAM_JUVEN=tany_rna_juveniles.bam

tar -xf /staging/lnell/ill/rna/${RNA_READS_JUVEN_TAR} -C ./
check_exit_status "extract rna juvenile" $?

gunzip ${RNA_READS_JUVEN1}.gz \
    && gunzip ${RNA_READS_JUVEN2}.gz
check_exit_status "gunzip rna - juveniles" $?

hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_JUVEN1} -2 ${RNA_READS_JUVEN2} \
    -k 3 -p ${THREADS} \
    | samtools sort -O bam - \
    > ${BAM_JUVEN}
check_exit_status "hisat2 - juveniles" $?

cp ${BAM_JUVEN} ${TARGET}/

rm ${RNA_READS_JUVEN1} ${RNA_READS_JUVEN2}

rm tany_hisat_idx*


conda deactivate





#' ===========================================================================
#' ===========================================================================
#'
#' Run BRAKER using RNAseq BAM files
#'
#' ===========================================================================
#' ===========================================================================


conda activate annotate-env


braker.pl --species=tanytarsus_gracilentus --genome=${GENOME} \
    --cores=${THREADS} \
    --softmasking \
    --bam=${BAM_ADULT},${BAM_JUVEN}

mv braker ${OUT_DIR}

# Saving output:
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working

