#!/bin/bash

#' Use RepeatModeler to create a repeat library for the Aedes aegypti assembly.
#'


time_stamp () {
    date +'%e %b %Y%t%T%n'
}
export progress_file=/staging/lnell/repeats-Aaegyp-progress.txt
echo -e "\n\n" > $progress_file
echo "-----------------------------------------------------------------------" \
    >> $progress_file
echo -e "\n\n" >> $progress_file
echo "           START RUN" >> $progress_file
time_stamp >> $progress_file
echo -e "\n\n" >> $progress_file
echo "-----------------------------------------------------------------------" \
    >> $progress_file
echo -e "\n\n" >> $progress_file



#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotations

export OUT_PREFIX=Aaegyp

export ASSEMBLY_FULL_PATH=/staging/lnell/assemblies/Aaegyp_contigs.fasta.gz


if [ ! -f "${ASSEMBLY_FULL_PATH}" ]; then
    echo "ERROR: Your assembly file ('${ASSEMBLY_FULL_PATH}') does not exist." 1>&2
    exit 1
fi
if ! [[ "${ASSEMBLY_FULL_PATH}" =~ (.fasta|.fa|.fasta.gz|.fa.gz)$ ]]; then
    echo -n "ERROR: Assembly must end in '.fasta', '.fa', '.fasta.gz', or '.fa.gz'. " 1>&2
    echo "Yours is '${ASSEMBLY_FULL_PATH}'." 1>&2
    exit 1
fi

if [ ! -d "${OUTPUT_LOC}" ]; then
    echo "ERROR: '${OUTPUT_LOC}' does not exist." 1>&2
    exit 1
fi





#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================



#' This should be divisible by 4 bc that's how many threads are used per
#' job using RMBlast, the default engine in RepeatModeler
export THREADS=$(count_threads)
if (( $THREADS < 4 )); then
    echo "ERROR: You need >= 4 threads, but you provided $THREADS." 1>&2
    exit 1
fi

. /app/.bashrc
conda activate repeat-env


export OUT_DIR=${OUT_PREFIX}_repeats_model
mkdir ${OUT_DIR}
cd ${OUT_DIR}


#' Naming outputs:
#' ------
#' Library of repeats from RepeatModeler:
export REPEATS_LIB=${OUT_PREFIX}_repeats_lib.fasta
#' Soft-masked assembly for use in BRAKER:
export MASKED_ASSEMBLY=${OUT_PREFIX}_contigs_masked.fasta
#' GFF of repeats in assembly:
export REPEAT_LOCS=${OUT_PREFIX}_repeats_locs.gff3


cp ${ASSEMBLY_FULL_PATH} ./
check_exit_status "cp genome" $?
export ASSEMBLY=$(basename "${ASSEMBLY_FULL_PATH}")
if [[ "${ASSEMBLY}" == *.gz ]]; then
    gunzip ${ASSEMBLY}
    check_exit_status "gunzip genome" $?
    ASSEMBLY=${ASSEMBLY%.gz}
fi



#' ------------------------------------------------
#' prep for RepeatModeler
#' ------------------------------------------------

conda activate main-env


mv ${ASSEMBLY} ${ASSEMBLY%.*}_orig.fasta
check_exit_status "rename genome" $?

#' Below does the following (in order):
#' - Simplify sequence names in headers (remove all after first space)
#' - Undo softmasking (I'll be manually softmasking repeat regions, if desired, later)
#' - Change to narrow (i.e., multi-line) format
sed '/^>/ s/ .*//' ${ASSEMBLY%.*}_orig.fasta \
    | tr [:lower:] [:upper:] \
    | seqtk seq -l 60 - \
    > ${ASSEMBLY}
rm ${ASSEMBLY%.*}_orig.fasta

conda deactivate


echo "finished RM prep" >> $progress_file
time_stamp >> $progress_file


#' ------------------------------------------------
#' run RepeatModeler
#' ------------------------------------------------

conda activate repeat-env


mkdir modeler
cd modeler
mv ../${ASSEMBLY} ./

# (RMBlast is default but put here for clarity)
BuildDatabase -engine rmblast -name ${OUT_PREFIX} ${ASSEMBLY} \
    >& BuildDatabase.log
check_exit_status "BuildDatabase" $?

echo "finished BuildDatabase" >> $progress_file
time_stamp >> $progress_file


echo "" > /staging/lnell/repeats-Aaegyp-RepeatModeler.out
echo "" > /staging/lnell/repeats-Aaegyp-RepeatModeler.err


timeout 65h \
    RepeatModeler -database ${OUT_PREFIX} -pa $(( THREADS / 4 )) -LTRStruct \
    -srand 633659157 \
    1> >(tee -a /staging/lnell/repeats-Aaegyp-RepeatModeler.out) \
    2> >(tee -a /staging/lnell/repeats-Aaegyp-RepeatModeler.err)
timeout_exit_status=$?
# check_exit_status "RepeatModeler" $?

cp /staging/lnell/repeats-Aaegyp-RepeatModeler.out RepeatModeler.out

cd ..


# Code 124 means it was still running
if [ $timeout_exit_status -eq 124 ]; then
    echo "ended RepeatModeler early" >> $progress_file
    time_stamp >> $progress_file
    cd ..
    tar -czf ${OUT_DIR}-unfinished.tar.gz ${OUT_DIR}
    mv ${OUT_DIR}-unfinished.tar.gz "${OUTPUT_LOC}"/
    rm -r ${OUT_DIR}
    # exit 85
    exit 1
fi



echo "finished RepeatModeler" >> $progress_file
time_stamp >> $progress_file





#' ------------------------------------------------
#' handling final output
#' ------------------------------------------------


# gzip < ${MASKED_ASSEMBLY} > ${MASKED_ASSEMBLY}.gz
# mv ${MASKED_ASSEMBLY}.gz "${OUTPUT_LOC}"/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz "${OUTPUT_LOC}"/
rm -r ${OUT_DIR}



exit 0
