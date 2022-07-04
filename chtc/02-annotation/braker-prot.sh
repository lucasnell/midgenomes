#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using proteins from a database (pipeline C)
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation

export OUT_DIR=tany_braker_prot

mkdir working
cd working

export GENOME=tany_contigs_masked.fasta


cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp genome" $?


# Using arthropoda proteins from OrthoDB v10:
wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar -xzf odb10_arthropoda_fasta.tar.gz
cat arthropoda/Rawdata/* > odb10_arthropoda.fasta
rm -r arthropoda odb10_arthropoda_fasta.tar.gz


conda activate annotate-env

braker.pl --genome=${GENOME} --prot_seq=odb10_arthropoda.fasta \
    --softmasking --cores=${THREADS}
check_exit_status "braker" $?

mv braker ${OUT_DIR}

# Saving output:
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working

