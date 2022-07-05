#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using RNAseq alignments (pipeline B)
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation

export OUT_DIR=Pstein_braker_rna

mkdir working
cd working

export GENOME=Pstein_contigs_masked.fasta

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


hisat2-build -q ${GENOME} Pstein_hisat_idx
check_exit_status "hisat2-build" $?

#' SRA accession numbers for these reads:
export SRA_NUMS=(SRR3951283 SRR3951284 SRR3951285)
#' Use above to construct array of BAM files:
export BAM_FILES=(${SRA_NUMS[@]/%/.bam})
export N_BAM_FILES=${#BAM_FILES[@]}

for ((i=0; i<${N_BAM_FILES}; i++)); do

    TAR_FILE=trimmed_Pstein_RNA_${SRA_NUMS[i]}.tar
    BAM_FILE=${BAM_FILES[i]}

    READS1=$(read_tar_name /staging/lnell/ill/rna/${TAR_FILE} 1)
    check_exit_status "extract read names - ${TAR_FILE} #1" $?
    READS2=$(read_tar_name /staging/lnell/ill/rna/${TAR_FILE} 2)
    check_exit_status "extract read names - ${TAR_FILE} #2" $?

    tar -xf /staging/lnell/ill/rna/${TAR_FILE} -C ./
    check_exit_status "extract rna - ${TAR_FILE}" $?

    gunzip ${READS1} && gunzip ${READS2}
    check_exit_status "gunzip rna - ${TAR_FILE}" $?

    READS1=${READS1/.gz/}
    READS2=${READS2/.gz/}

    hisat2 -x Pstein_hisat_idx -1 ${READS1} -2 ${READS2} -k 3 -p ${THREADS} \
        | samtools sort -O bam - \
        > ${BAM_FILE}
    check_exit_status "hisat2 - adults" $?

    cp ${BAM_FILE} ${TARGET}/

    rm ${READS1} ${READS2}
    unset READS1 READS2 TAR_FILE BAM_FILE

done


rm Pstein_hisat_idx*

conda deactivate





#' ===========================================================================
#' ===========================================================================
#'
#' Run BRAKER using RNAseq BAM files
#'
#' ===========================================================================
#' ===========================================================================


conda activate annotate-env


braker.pl --species=Parochlus_steinenii --genome=${GENOME} \
    --cores=${THREADS} \
    --softmasking \
    --bam=$(IFS=, ; echo "${BAM_FILES[*]}")

mv braker ${OUT_DIR}

# Saving output:
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working

