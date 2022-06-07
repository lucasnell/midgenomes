#!/bin/bash

#'
#' Polish assembly using ONT reads with Racon and medaka.
#' Requires argument for which assembly to polish.
#'

export ASSEMBLY=$1

if [ -z ${2+x} ]; then
    N_RACON=3
else
    N_RACON=$2
    if ! [[ $N_RACON =~ ^[0-9]+$ ]]; then
        echo "ERROR: Second arg is not an integer! " 1>&2
        exit 1
    fi
fi
export N_RACON



if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi


. /app/.bashrc
conda activate assembly-env

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
if (( N_RACON > 0 )); then
    OUT_NAME=${ASSEMBLY/.fasta/}_ont-polish
else
    OUT_NAME=${ASSEMBLY/.fasta/}_medaka
fi
export OUT_NAME
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

export LONGREADS=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?


# ===================*
# Racon
# ===================*


export ALIGN=mm2_align.paf
export NEW_ASSEMBLY=${ASSEMBLY}

# X rounds of polishing using Racon:
if (( N_RACON > 0 )); then
    for i in $(seq 1 ${N_RACON}); do
        minimap2 -x map-ont -t $((THREADS - 1)) -K 1G \
            ${NEW_ASSEMBLY} ${LONGREADS} \
            > ${ALIGN}
        # These parameters below are for compatibility with medaka:
        # -m 8 -x -6 -g -8 -w 500
        racon -m 8 -x -6 -g -8 -w 500 -t ${THREADS} \
            ${LONGREADS} ${ALIGN} ${NEW_ASSEMBLY} \
            > racon_${i}.fasta
        NEW_ASSEMBLY=racon_${i}.fasta
    done
    rm ${ALIGN}

    # Summarize consensus after last round:

    summ-scaffs.py ${NEW_ASSEMBLY} | tee racon_contigs_summary.out
    check_exit_status "summ-scaffs.py (racon)" $?

    run_busco ${NEW_ASSEMBLY} ${THREADS}
    mv busco.out racon_busco.out
    rm -r busco_downloads busco

    pretty-csv.py -s racon_contigs_summary.out -b racon_busco.out \
        ${OUT_NAME}_racon \
        | tee ${OUT_NAME}_racon.csv

fi




# ===================*
# Medaka
# ===================*

#' For the model to use (`-m` argument), the "r941_min_" should remain
#' constant even if you use a later version of guppy because this specifies
#' the pore and device.
medaka_consensus -i ${LONGREADS} -d ${NEW_ASSEMBLY} -o medaka -t ${THREADS} \
    -m r941_min_hac_g507

cp ./medaka/consensus.fasta ${OUT_FASTA}

# Summarize medaka consensus:

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads busco


pretty-csv.py -s contigs_summary.out -b busco.out ${OUT_NAME} \
    | tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

rm ${LONGREADS} ${ASSEMBLY}

# These files are massive. Removing for storage.
cd medaka
calls_to_draft.bam
rm calls_to_draft.* consensus_probs.hdf
cd ..

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

