#!/bin/bash

#'
#' Polish assembly using Illumina reads with NextPolish.
#' Requires argument for which assembly to polish.
#'

#' Not immediately exiting so that all errors are printed.
export status=0

#' Requires 3 arguments:
if [[ $# -ne 3 ]]; then
    echo "ERROR: nextpolish.sh requires 3 input arguments. You have $#." 1>&2
    status=1
fi

#' First is assembly to polish:
export ASSEMBLY=$1
if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo -n "ERROR: First arg to nextpolish.sh must end in *.fasta. " 1>&2
    echo "Yours is '${ASSEMBLY}'." 1>&2
    status=1
fi

#' Second argument is number of rounds.
#' I wouldn't recommend more than 3: More rounds isn't necessarily better!
export N_ROUNDS=$2
if ! [[ $N_ROUNDS =~ ^[0-9]+$ ]]; then
    echo -n "ERROR: Second arg to nextpolish.sh must be an integer. " 1>&2
    echo "Your is '${N_ROUNDS}'." 1>&2
    status=1
fi

#' Third argument is whether Illumina reads should be normalized ($NORM == 1)
#' or not ($NORM == 0).
export NORM=$3
if ! [[ "${NORM}" == "0" ]] && ! [[ "${NORM}" == "1" ]]; then
    echo "ERROR: Third arg to nextpolish.sh should be 0 or 1. " 1>&2
    echo "Yours is '${NORM}'." 1>&2
    status=1
fi

if (( status != 0 )); then
    exit 1
else
    echo -e "\nBasic tests passed for nextpolish.sh inputs.\n\n"
fi


. /app/.bashrc
conda activate main-env
source /staging/lnell/helpers.sh


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
if (( NORM == 0 )); then
    OUT_NAME=${ASSEMBLY/.fasta/}_nextpolish
else
    OUT_NAME=${ASSEMBLY/.fasta/}_norm-nextpolish
fi
export OUT_NAME
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

#' Input assembly:
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

#' All reads have uncalled bases filtered out (to prevent N in final assembly),
#' and if NORM == 1, then they are also normalized.

READS1=noN_trimmed_MyKS-19-B_S18_L002_R1_001.fastq
READS2=noN_trimmed_MyKS-19-B_S18_L002_R2_001.fastq
READS_TAR=noN_trimmed_MyKS-19-B_S18.tar.gz
if (( NORM == 1 )); then
    READS1=norm_${READS1}
    READS2=norm_${READS2}
    READS_TAR=norm_${READS_TAR}
fi
export READS1
export READS2
export READS_TAR

cp /staging/lnell/dna/trimmed/for-polishing/${READS_TAR} ./ \
    && tar -xzf ${READS_TAR} \
    && rm ${READS_TAR}


#' For proper date-times:
export TZ="America/Chicago"

#' Function to print to stdout and stderr.
#' Usage:
#'     print_start [NAME_OF_STAGE]
print_start () {
    local L="================================================================\n"
    echo -e "\n\n\n${L}${L}\n""$@" "\n@" $(date "+%F %T") "\n\n${L}${L}\n" \
        | tee /dev/stderr
}
#' Similar but for end point. No arguments needed here.
print_end () {
    local L="================================================================\n"
    echo -e "\n\n\n\n${L}${L}\nFinished\n@ "$(date "+%F %T") "\n\n" \
        | tee /dev/stderr
}

#'
#' Functions for each section of NextPolish pipeline, to easy capturing output.
#' All these functions use environmental variables, not input arguments.
#'

status=0
#' Index assembly then map reads to it.
do_map () {
    print_start "BWA-indexing assembly"
    bwa index ${input}
    status=$?
    if (( status != 0 )); then return $status; fi
    print_start "Aligning"
    bwa mem -t ${THREADS} ${input} ${READS1} ${READS2} | \
        samtools view --threads 3 -F 0x4 -b -| \
        samtools fixmate -m --threads 3 - -| \
        samtools sort -m 2g --threads 5 -| \
        samtools markdup --threads 5 -r - sgs.sort.bam
    status=$?
    print_end
    return $status
}
#' Index BAM and faidx-index assembly.
do_index () {
    # print_start "Index BAM"
    samtools index -@ ${THREADS} sgs.sort.bam
    status=$?
    if (( status != 0 )); then return $status; fi
    # print_start "FAIDX-index assembly"
    samtools faidx ${input}
    status=$?
    # print_end
    return $status
}
#' Polish assembly using NextPolish, remove unneeded files from previous steps,
#' then update `input` variable.
#' NOTE: If you change the iterator for step from `j` to something else,
#'       you should also update this function!
do_polish () {
    print_start "Polishing"
    local OUT=genome_round${i}_step${j}.fa
    python /opt/NextPolish/lib/nextpolish1.py -g ${input} \
        --task ${j} -p ${THREADS} -s sgs.sort.bam \
        > ${OUT}
    status=$?
    if (( status != 0 )); then return $status; fi
    rm ${input}.* sgs.sort.bam*
    input=${OUT}
    print_end
    return $status
}


export input=${ASSEMBLY}

mkdir logs

for ((i=1; i<=${N_ROUNDS}; i++)); do
    for ((j=1; j<=2; j++)); do
        # Map:
        do_map 2> logs/mapping_round${i}_step${j}.log
        # Only print this massive *.log file if something went wrong:
        if (( status != 0 )); then
            cat logs/mapping_round${i}_step${j}.log >&2
        fi
        check_exit_status "alignment, round ${i}, step ${j}" $status
        # Index:
        do_index
        check_exit_status "index, round ${i}, step ${j}" $status
        # Polish:
        do_polish 2> >(tee -a logs/polish_round${i}_step${j}.log >&2)
        check_exit_status "polish, round ${i}, step ${j}" $status
    done
done


# Finally polished genome file: $input
cp ${input} ${OUT_FASTA}
check_exit_status "cp output" $?

rm ${READS1} ${READS2} ${ASSEMBLY}

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
