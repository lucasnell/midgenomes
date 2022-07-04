#!/bin/bash

# Use RepeatModeler to create a repeat library for the Tanytarsus gracilentus genome
# and to annotate the genome based on this library.


#' This should be divisible by 4 bc that's how many threads are used per
#' job using RMBlast, the default engine in RepeatModeler
export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


. /app/.bashrc
conda activate repeat-env


export OUT_DIR=tany_repeats

mkdir ${OUT_DIR}
cd ${OUT_DIR}


# Soft-masked assembly for use in BRAKER:
export MASKED_ASSEMBLY=tany_contigs_masked.fasta
# GFF of repeats in assembly:
export REPEAT_LOCS=tany_repeats_locs.gff3
# Library of repeats from RepeatModeler:
export REPEATS_LIB=tany_repeats_lib.fasta

export GENOME=tany_contigs.fasta
cp /staging/lnell/assemblies/${GENOME}.gz ./ && gunzip ${GENOME}.gz




#' ------------------------------------------------
#' run RepeatModeler
#' ------------------------------------------------

mkdir modeler
cd modeler
mv ../${GENOME} ./

# (RMBlast is default but put here for clarity)
BuildDatabase -engine rmblast -name tgrac ${GENOME} \
    >& BuildDatabase.log
check_exit_status "BuildDatabase" $?

RepeatModeler -database tgrac -pa $(( THREADS / 4 )) -LTRStruct \
    1> >(tee -a RepeatModeler.out)
check_exit_status "RepeatModeler" $?





#' ------------------------------------------------
#' softmask assembly using RepeatMasker
#' ------------------------------------------------

mkdir masker
cd masker

mv ../modeler/${GENOME} ./
cp ../modeler/tgrac-families.fa ${REPEATS_LIB}

RepeatMasker \
    -s -xsmall -gff \
    -lib ${REPEATS_LIB} \
    -pa $(( THREADS / 4 )) \
    ${GENOME} \
    1> >(tee -a RepeatMasker.out)


cp ${GENOME}.masked ../${MASKED_ASSEMBLY}
cp ${GENOME}.out.gff ../${REPEAT_LOCS}
mv ${REPEATS_LIB} ../
rm ${GENOME}
cd ..



#' ------------------------------------------------
#' handling output
#' ------------------------------------------------


gzip < ${MASKED_ASSEMBLY} > ${MASKED_ASSEMBLY}.gz


mv ${MASKED_ASSEMBLY}.gz /staging/lnell/annotation/


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/annotation/
rm -r ${OUT_DIR}


