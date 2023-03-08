#!/bin/bash

#' Use RepeatModeler to create a repeat library for the Aedes aegypti assembly.
#' This recovers from a run that had to be ended early.
#'
#' Takes a single input, which is the full path to a *.tar.gz file output from
#' repeats-Aaegyp-model.sh
#'


time_stamp () {
    date +'%e %b %Y%t%T%n'
}

export start_date=$(date +'%e-%b-%Y')
export progress_file=/staging/lnell/repeats-Aaegyp-recover-${start_date}-progress.txt
export std_file_prefix=repeats-Aaegyp-recover-${start_date}
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


export INPUT_TAR=$1
if [ ! -f "${INPUT_TAR}" ]; then
    echo "ERROR: Your input tar file ('${INPUT_TAR}') does not exist." 1>&2
    exit 1
fi

export OUTPUT_LOC=/staging/lnell/annotations
export OUT_PREFIX=Aaegyp
export ASSEMBLY=Aaegyp_contigs.fasta

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
export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')
if (( $THREADS < 4 )); then
    echo "ERROR: You need >= 4 threads, but you provided $THREADS." 1>&2
    exit 1
fi

. /app/.bashrc
conda activate repeat-env


export OUT_DIR=${OUT_PREFIX}_repeats_model




#' Naming outputs:
#' ------
#' Library of repeats from RepeatModeler:
export REPEATS_LIB=${OUT_PREFIX}_repeats_lib.fasta
#' Soft-masked assembly for use in BRAKER:
export MASKED_ASSEMBLY=${OUT_PREFIX}_contigs_masked.fasta
#' GFF of repeats in assembly:
export REPEAT_LOCS=${OUT_PREFIX}_repeats_locs.gff3




#' ------------------------------------------------
#' Output from previous RepeatModeler step
#' ------------------------------------------------


tar -xf ${INPUT_TAR} -C ./

if [ ! -d "${OUT_DIR}" ]; then
    echo "ERROR: '${INPUT_TAR}' is not un-tarring to ${OUT_DIR}." 1>&2
    exit 1
fi
cd ${OUT_DIR}
if [ ! -d "modeler" ]; then
    echo "ERROR: '${INPUT_TAR}' does not contain 'modeler' folder." 1>&2
    exit 1
fi




#' ------------------------------------------------
#' run RepeatModeler
#' ------------------------------------------------

conda activate repeat-env


cd modeler

if [ ! -f "${ASSEMBLY}" ]; then
    echo "ERROR: No assembly file ('${ASSEMBLY}') inside 'modeler' folder." 1>&2
    exit 1
fi
# It's assumed that BuildDatabase was already run
if [ ! -f "${OUT_PREFIX}.nhr" ]; then
    echo "ERROR: No BuildDatabase file '${OUT_PREFIX}.nhr' inside 'modeler' folder." 1>&2
    exit 1
fi

# Find recovery folder:
if (( $(find . -type d -name "RM_*") > 1 )); then
    echo "ERROR: Multiple recovery folders found." 1>&2
    exit 1
else
    export recover_dir=$(find . -type d -name "RM_*" | sed 's|^\./||g')
fi

#' To allow recovery of RepeatModeler after round 6, but before LTRStruct
#' See https://github.com/Dfam-consortium/RepeatModeler/issues/80#issuecomment-617265216
if [ -f "${recover_dir}/round-6/consensi.fa" ]; then
    mv ${recover_dir}/round-6/consensi.fa ${recover_dir}/round-6/consensi.fa.bak
fi

timeout 65h \
    RepeatModeler -database ${OUT_PREFIX} -pa $(( THREADS / 4 )) -LTRStruct \
    -srand 633659157 -recoverDir ${recover_dir} \
    1> >(tee -a /staging/lnell/${std_file_prefix}.out) \
    2> >(tee -a /staging/lnell/${std_file_prefix}.err)
timeout_exit_status=$?
# check_exit_status "RepeatModeler" $?

cp /staging/lnell/${std_file_prefix}.out RepeatModeler.out

cd ..


#' Code 124 means it was still running when time limit cut it off
if [ $timeout_exit_status -eq 124 ]; then
    echo "ended RepeatModeler early" >> $progress_file
    time_stamp >> $progress_file
    cd ..
    tar -czf ${OUT_DIR}-unfinished-${start_date}.tar.gz ${OUT_DIR}
    mv ${OUT_DIR}-unfinished-${start_date}.tar.gz "${OUTPUT_LOC}"/
    rm -r ${OUT_DIR}
    # exit 85
    exit 1
fi



# LEFT OFF HERE



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



#' If this continues to fail, see these issues:
#' https://github.com/Dfam-consortium/RepeatModeler/issues/80
#' https://github.com/Dfam-consortium/RepeatModeler/issues/69

exit 0
