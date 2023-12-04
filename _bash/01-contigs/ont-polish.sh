#!/bin/bash

#'
#' Polish assembly using ONT reads with Racon and medaka.
#' Requires argument for which assembly to polish.
#' Note: options MUST come before assembly.
#'
#' Usage:
#' ont-polish.sh [options] ASSEMBLY
#'
#' Options:
#'   -n Nanopore reads to use. Must be in `/staging/lnell/ont` and end
#'      with `.fastq`. Defaults to `basecalls_guppy-5.0.11.fastq`.
#'   -r Number of rounds of Racon polishing. Must be an integer >= 0.
#'      Defaults to 3.
#'   -m Number of rounds of medaka polishing. Must be 0 or 1. Defaults to 1.
#'


export LONGREADS=basecalls_guppy-5.0.11.fastq
export N_RACON=3
export N_MEDAKA=1


while getopts ":n:r:m:" opt; do
    case $opt in
        n)
            LONGREADS="$OPTARG"
            if [[ "${LONGREADS}" != *.fastq ]]; then
                echo -n "ERROR: -n arg must end in *.fastq. " 1>&2
                echo "Yours is '${LONGREADS}'." 1>&2
                exit 1
            fi
            ;;
        r)
            N_RACON="$OPTARG"
            if ! [[ $N_RACON =~ ^[0-9]+$ ]]; then
                echo -n "ERROR: -r arg should be an integer. " 1>&2
                echo "Yours is '${N_RACON}'." 1>&2
                exit 1
            fi
            ;;
        m)
            N_MEDAKA="$OPTARG"
            if [[ $N_MEDAKA != "0" ]] && [[ $N_MEDAKA != "1" ]]; then
                echo -n "ERROR: -m arg should be 0 or 1. " 1>&2
                echo "Yours is '${N_MEDAKA}'." 1>&2
                exit 1
            fi
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            exit 1
            ;;
    esac
done

export ASSEMBLY=${@:$OPTIND:1}
if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo -n "ERROR: Assembly must end in *.fasta. " 1>&2
    echo "Yours is '${ASSEMBLY}'." 1>&2
    exit 1
fi


if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi

if (( N_RACON + N_MEDAKA == 0 )); then
    echo "Both Racon and medaka cannot have zero rounds." 1>&2
    exit 1
fi



. /app/.bashrc
conda activate assembly-env

export THREADS=$(count_threads)

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
if (( N_RACON > 0 )) && ((N_MEDAKA > 0)); then
    OUT_NAME=${ASSEMBLY/.fasta/}_ont-polish
elif (( N_RACON > 0 )); then
    OUT_NAME=${ASSEMBLY/.fasta/}_racon
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

cp /staging/lnell/ont/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp long reads" $?


# ===================*
# Racon
# ===================*


export NEW_ASSEMBLY=${ASSEMBLY}

# X rounds of polishing using Racon:
if (( N_RACON > 0 )); then
    ALIGN=mm2_align.paf
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

    #' Summarize consensus after last round, but if no medaka will be run, then
    #' we can just do this at the end.
    if (( N_MEDAKA > 0 )); then

        summ-scaffs.py ${NEW_ASSEMBLY} | tee racon_contigs_summary.out
        check_exit_status "summ-scaffs.py (racon)" $?

        run_busco ${NEW_ASSEMBLY} ${THREADS}
        mv busco.out racon_busco.out
        rm -r busco_downloads busco

        pretty-csv.py -s racon_contigs_summary.out -b racon_busco.out \
            ${OUT_NAME}_racon \
            | tee ${OUT_NAME}_racon.csv

    else

        cp racon_${i}.fasta ${OUT_FASTA}

    fi

fi




# ===================*
# Medaka
# ===================*

if (( N_MEDAKA > 0 )); then

    #' (The note below is for Tanytarsus gracilentus assembly only)
    #' For the model to use (`-m` argument), the "r941_min_" should remain
    #' constant even if you use a later version of guppy because this specifies
    #' the pore and device.
    medaka_consensus -i ${LONGREADS} -d ${NEW_ASSEMBLY} -o medaka -t ${THREADS} \
        -m r941_min_hac_g507

    cp ./medaka/consensus.fasta ${OUT_FASTA}

    # These files are massive. Removing for storage.
    cd medaka
    rm calls_to_draft.* consensus_probs.hdf
    cd ..

fi


# Summarize final consensus:

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

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

