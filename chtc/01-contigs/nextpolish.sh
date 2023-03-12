#!/bin/bash


#'
#' Polish assembly using Illumina reads with NextPolish.
#' Requires argument for which assembly to polish.
#' Note: options MUST come before assembly.
#'
#' Usage:
#' nextpolish.sh [options] ASSEMBLY
#'
#' Options:
#'   -i tar file containing Illumina reads to use.
#'      Must be in `/staging/lnell/ill/dna/trimmed` and end with `.tar`.
#'      This script also assumes that only reads are inside the tar file and
#'      that if the file names are sorted alphabetically, the first file is
#'      #1 of pair, and the second is #2 of pair.
#'      Defaults to `trimmed_MyKS-19-B_S18.tar`.
#'   -r Number of rounds of polishing. Must be an integer > 0.
#'      Defaults to 3.
#'

export READS_TAR=trimmed_MyKS-19-B_S18.tar
export N_ROUNDS=3
export READS_LOC=/staging/lnell/ill/dna/trimmed

while getopts ":i:r:" opt; do
    case $opt in
        i)
            READS_TAR="$OPTARG"
            if [[ "${READS_TAR}" != *.tar ]]; then
                echo -n "ERROR: -n arg must end in *.tar. " 1>&2
                echo "Yours is '${READS_TAR}'." 1>&2
                exit 1
            fi
            if [ ! -f ${READS_LOC}/${READS_TAR} ]; then
                echo -n "ERROR: '${READS_LOC}/${READS_TAR}' does not exist." 1>&2
                exit 1
            fi
            ;;
        r)
            N_ROUNDS="$OPTARG"
            if ! [[ $N_ROUNDS =~ ^[0-9]+$ ]] || (( N_ROUNDS <= 0 )); then
                echo -n "ERROR: -r arg should be an integer > 0. " 1>&2
                echo "Yours is '${N_ROUNDS}'." 1>&2
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


#' Extract individual read names from tar file:
export READS1=$(read_tar_name ${READS_LOC}/${READS_TAR} 1)
export READS2=$(read_tar_name ${READS_LOC}/${READS_TAR} 2)



. /app/.bashrc
conda activate main-env


export THREADS=$(count_threads)


# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_nextpolish
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

#' Input assembly:
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

#' Function to print message with a date-time afterward.
#' Used below to make the log files more readable.
#' Usage:
#'     msg [NAME_OF_STAGE]
msg () {
    echo -e "\n\n$@\n@" $(date "+%F %T") "\n\n"
}




#' ===========================================================================
#' ===========================================================================
#'
#' Filter uncalled bases from reads to prevent N in final assembly:
#'
#' ===========================================================================
#' ===========================================================================

# Input reads that have already been trimmed for adapters, etc.
tar -xf ${READS_LOC}/${READS_TAR} -C ./
check_exit_status "extract reads tar file" $?

export NON_READS1=noN_${READS1%.gz}
export NON_READS2=noN_${READS2%.gz}

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${NON_READS1} --out2 ${NON_READS2} \
    --thread ${THREADS} \
    --disable_length_filtering \
    --disable_adapter_trimming \
    --disable_trim_poly_g \
    --dont_eval_duplication \
    --qualified_quality_phred 0 \
    --unqualified_percent_limit 100 \
    --n_base_limit 0
rm ${READS1} ${READS2} fastp*

READS1=${NON_READS1}
READS2=${NON_READS2}

unset NON_READS1 NON_READS2 READS_TAR




#' ===========================================================================
#' ===========================================================================
#'
#' Polish using NextPolish
#'
#' ===========================================================================
#' ===========================================================================

#' This variable will be changed for each iteration below, but it obviously
#' starts as the input assembly.
export CURRENT_ASS=${ASSEMBLY}

#' Temporary bam file
export S_BAM=sorted.bam


for ((i=1; i<=${N_ROUNDS}; i++)); do
    for ((j=1; j<=2; j++)); do

        # ---------------------
        # Map:
        LOG=mapping_round${i}_step${j}.log
        msg "BWA-indexing assembly" > ${LOG}
        bwa index ${CURRENT_ASS} 2>> ${LOG}
        check_exit_status "bwa index, round ${i}, step ${j}" $?
        msg "Aligning" >> ${LOG}
        bwa mem -t ${THREADS} ${CURRENT_ASS} ${READS1} ${READS2} \
            | samtools view --threads 3 -F 0x4 -b - \
            | samtools fixmate -m --threads 3 - - \
            | samtools sort -m 2g --threads 5 - \
            | samtools markdup --threads 5 -r - ${S_BAM} \
            2>> ${LOG}
        status=$?
        # Only print this massive *.log file if something went wrong:
        if (( status != 0 )); then cat ${LOG} >&2; fi
        check_exit_status "alignment, round ${i}, step ${j}" $status
        msg "Finished" >> ${LOG}

        # ---------------------
        # Index:
        samtools index -@ ${THREADS} ${S_BAM}
        check_exit_status "samtools index, round ${i}, step ${j}" $?
        samtools faidx ${CURRENT_ASS}
        check_exit_status "samtools faidx, round ${i}, step ${j}" $?

        # ---------------------
        # Polish:
        LOG=polish_round${i}_step${j}.log
        OUT=genome_round${i}_step${j}.fa
        msg "Polishing" > ${LOG}
        python /opt/NextPolish/lib/nextpolish1.py -g ${CURRENT_ASS} \
            --task ${j} -p ${THREADS} -s ${S_BAM} \
            > ${OUT} \
            2>> >(tee -a ${LOG} >&2)
        check_exit_status "polish, round ${i}, step ${j}" $?
        rm ${CURRENT_ASS}.* ${S_BAM}*
        msg "Finished" >> ${LOG}
        CURRENT_ASS=${OUT}

        unset LOG OUT

    done
done




#' ===========================================================================
#' ===========================================================================
#'
#' Handle output
#'
#' ===========================================================================
#' ===========================================================================

cp ${CURRENT_ASS} ${OUT_FASTA}
check_exit_status "cp output" $?

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco busco_downloads

pretty-csv.py -s contigs_summary.out -b busco.out ${OUT_NAME} \
    | tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

rm ${READS1} ${READS2} ${ASSEMBLY}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


