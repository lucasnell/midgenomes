#!/bin/bash

#'
#' Polish assembly using Nanopore and Illumina reads with NextPolish.
#' Requires argument for which assembly to polish and number of rounds for
#' long and short reads.
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

#' Second argument is number of rounds for long reads.
export LONG_ROUNDS=$2
if ! [[ $LONG_ROUNDS =~ ^[0-9]+$ ]]; then
    echo -n "ERROR: Second arg to nextpolish.sh must be an integer >= 0. " 1>&2
    echo "Your is '${LONG_ROUNDS}'." 1>&2
    status=1
fi

#' Third argument is number of rounds for short reads.
#' I wouldn't recommend more than 3: More rounds isn't necessarily better!
export SHORT_ROUNDS=$2
if ! [[ $SHORT_ROUNDS =~ ^[0-9]+$ ]]; then
    echo -n "ERROR: Third arg to nextpolish.sh must be an integer >= 0. " 1>&2
    echo "Your is '${SHORT_ROUNDS}'." 1>&2
    status=1
fi

#' If any of above failed, exit
if (( status != 0 )); then
    exit 1
else
    echo -e "\nBasic tests passed for nextpolish.sh inputs.\n\n"
fi

#' Can't have zeros for both numbers of rounds.
if (( LONG_ROUNDS == 0 )) && (( SHORT_ROUNDS == 0 )); then
    echo "ERROR: Second or third arg to nextpolish.sh must be > 0." 1>&2
    exit 1
fi


# export ASSEMBLY=contigs_quickmerge_medaka.fasta
# export LONG_ROUNDS=3
# export SHORT_ROUNDS=3


. /app/.bashrc
conda activate main-env
source /staging/lnell/helpers.sh


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_nextpolish-${LONG_ROUNDS}-${SHORT_ROUNDS}
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
    echo -e "\n\n$@\n@" $(TZ=America/Chicago date "+%F %T") "\n\n"
}

#' This variable will be changed for each iteration below, but it obviously
#' starts as the input assembly.
export CURRENT_ASS=${ASSEMBLY}



#' ============================================================================
#' ============================================================================
#'
#' Polish with Nanopore reads:
#'
#' ============================================================================
#' ============================================================================




if (( LONG_ROUNDS > 0 )); then

    export RAW_LONG_READS=basecalls_guppy-5.0.11.fastq
    cp /staging/lnell/${RAW_LONG_READS}.gz ./ && gunzip ${RAW_LONG_READS}.gz
    check_exit_status "cp long reads" $?

    #' This hangs if I use all the data, so I filter by read length
    #' to reduce overall coverage.
    #'
    #' These are the thresholds for a number of coverages:
    #'   coverage thresh
    #'        100   9038
    #'         50  12202
    #'         40  12999
    #'         30  13963
    #'         20  15244
    #'         10  17282
    #'
    #' Code to create this table, where "sequencing_summary.txt.gz" is the
    #' tab-delimited summary table output from the ONT sequencer.
    #'
    #' ```r
    #' library(readr)
    #' lens <- read_tsv(paste0("sequencing_summary.txt.gz"))
    #' lens <- lens[["sequence_length_template"]]
    #' get_cov_diff <- function(thresh, desired, penalty = 1) {
    #'     observed <- sum(lens[lens > thresh]) / 92e6L
    #'     cov_diff <- abs(observed - desired)
    #'     # Add small penalty for being below desired:
    #'     if (observed < desired) cov_diff <- cov_diff + penalty
    #'     return(cov_diff)
    #' }
    #' thresh_df <- data.frame(coverage = c(100, 5:1 * 10),
    #'                         thresh = 0.0)
    #' set.seed(1703704117)
    #' thresh_df$thresh <- sapply(thresh_df$coverage,
    #'                            function(z) {
    #'                                fx <- function(x) get_cov_diff(x, z)
    #'                                opt = optimize(fx, c(1, 20e3L))
    #'                                return(floor(opt[["minimum"]]))
    #'                            })
    #' print(thresh_df, row.names = FALSE)
    #' ```
    export LONG_READS=filtered_nanopore_reads.fastq
    reformat.sh in=${RAW_LONG_READS} out=${LONG_READS} qin=33 minlength=9038

    # rm ${RAW_LONG_READS}

    export L_BAM=long_sorted.bam

    mkdir long_logs

    for ((i=1; i<=${LONG_ROUNDS}; i++)); do

        # ---------------------
        # Map:
        LOG=long_logs/mapping_round${i}.log
        msg "Aligning" > ${LOG}
        minimap2 -ax map-ont -t ${THREADS} ${CURRENT_ASS} ${LONG_READS} \
            | samtools sort - -m 2g --threads 5 -o ${L_BAM} \
            2>> >(tee -a ${LOG} >&2)
        check_exit_status "alignment, longs round ${i}" $?
        msg "Finished" >> ${LOG}

        # ---------------------
        # Index
        samtools index ${L_BAM}
        check_exit_status "samtools index, longs round ${i}" $?
        echo $(pwd)/${L_BAM} > ${L_BAM}.fofn

        # ---------------------
        # Polish
        #' (Have to go back to base environment bc `nextpolish2.py` requires
        #'  the `psutil` package.)
        conda deactivate
        LOG=long_logs/polish_round${i}.log
        OUT=genome_longs_round${i}.fa
        msg "Polishing" > ${LOG}
        python /opt/NextPolish/lib/nextpolish2.py -g ${CURRENT_ASS} \
            -l ${L_BAM}.fofn -r ont -p ${THREADS} -sp -a \
            > ${OUT} \
            2>> >(tee -a ${LOG} >&2)
        check_exit_status "polish, longs round ${i}" $?
        rm ${L_BAM}*
        msg "Finished" >> ${LOG}
        CURRENT_ASS=${OUT}
        conda activate main-env

        unset LOG OUT

    done

    rm ${LONG_READS}

    unset RAW_LONG_READS LONG_READS L_BAM

fi




#' ============================================================================
#' ============================================================================
#'
#' Polish with Illumina reads:
#'
#' ============================================================================
#' ============================================================================




if (( SHORT_ROUNDS > 0 )); then

    #' These Illumina reads have uncalled bases filtered out to prevent
    #' N in final assembly.
    export SHORT_READS1=noN_trimmed_MyKS-19-B_S18_L002_R1_001.fastq
    export SHORT_READS2=noN_trimmed_MyKS-19-B_S18_L002_R2_001.fastq
    export SHORT_READS_TAR=noN_trimmed_MyKS-19-B_S18.tar.gz

    cp /staging/lnell/dna/trimmed/for-polishing/${SHORT_READS_TAR} ./ \
        && tar -xzf ${SHORT_READS_TAR} \
        && rm ${SHORT_READS_TAR}
    check_exit_status "cp short reads" $?

    export S_BAM=short_sorted.bam

    mkdir short_logs

    for ((i=1; i<=${SHORT_ROUNDS}; i++)); do
        for ((j=1; j<=2; j++)); do

            # ---------------------
            # Map:
            LOG=short_logs/mapping_round${i}_step${j}.log
            msg "BWA-indexing assembly" > ${LOG}
            bwa index ${CURRENT_ASS} 2>> ${LOG}
            check_exit_status "bwa index, shorts round ${i}, step ${j}" $?
            msg "Aligning" >> ${LOG}
            bwa mem -t ${THREADS} ${CURRENT_ASS} ${SHORT_READS1} ${SHORT_READS2} \
                | samtools view --threads 3 -F 0x4 -b - \
                | samtools fixmate -m --threads 3 - - \
                | samtools sort -m 2g --threads 5 - \
                | samtools markdup --threads 5 -r - ${S_BAM} \
                2>> ${LOG}
            status=$?
            # Only print this massive *.log file if something went wrong:
            if (( status != 0 )); then cat ${LOG} >&2; fi
            check_exit_status "alignment, shorts round ${i}, step ${j}" $status
            msg "Finished" >> ${LOG}


            # ---------------------
            # Index:
            samtools index -@ ${THREADS} ${S_BAM}
            check_exit_status "samtools index, shorts round ${i}, step ${j}" $?
            samtools faidx ${CURRENT_ASS}
            check_exit_status "samtools faidx, shorts round ${i}, step ${j}" $?

            # ---------------------
            # Polish:
            LOG=short_logs/polish_round${i}_step${j}.log
            OUT=genome_round${i}_step${j}.fa
            msg "Polishing" > ${LOG}
            python /opt/NextPolish/lib/nextpolish1.py -g ${CURRENT_ASS} \
                --task ${j} -p ${THREADS} -s ${S_BAM} \
                > ${OUT} \
                2>> >(tee -a ${LOG} >&2)
            check_exit_status "polish, shorts round ${i}, step ${j}" $?
            rm ${CURRENT_ASS}.* ${S_BAM}*
            msg "Finished" >> ${LOG}
            CURRENT_ASS=${OUT}

            unset LOG OUT

        done
    done

    unset S_BAM SHORT_READS1 SHORT_READS2 SHORT_READS_TAR


fi



#' ============================================================================
#' ============================================================================
#'
#' Summarize, move final output:
#'
#' ============================================================================
#' ============================================================================



# Finally polished genome file: $CURRENT_ASS
cp ${CURRENT_ASS} ${OUT_FASTA}
check_exit_status "cp output" $?


















# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# # THIS VERSION SEEMS TO HANG ON THE CLUSTER...
#
# conda deactivate
#
# ls ${READS1} ${READS2} > sgs.fofn
# ls ${LONGREADS} > lgs.fofn
#
# cat << EOF > run.cfg
# [General]
# job_type = local
# job_prefix = nextPolish
# task = best
# rewrite = yes
# rerun = 3
# parallel_jobs = $(( THREADS / 8 ))
# multithread_jobs = 8
# genome = ./${ASSEMBLY}
# genome_size = auto
# workdir = ./01_rundir
# polish_options = -p {multithread_jobs}
#
# [sgs_option]
# sgs_fofn = ./sgs.fofn
# sgs_options = -max_depth 100 -bwa
#
# [lgs_option]
# lgs_fofn = ./lgs.fofn
# lgs_options = -min_read_len 1k -max_depth 100
# lgs_minimap2_options = -x map-ont -t {multithread_jobs}
# EOF
#
#
# nextPolish run.cfg
#
#
# cp ./01_rundir/genome.nextpolish.fasta ${OUT_FASTA}
#
# # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




# check_exit_status "cp output" $?
#
# rm ${READS1} ${READS2} ${LONGREADS} ${ASSEMBLY}
#
# conda activate main-env
#
# summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
# check_exit_status "summ-scaffs.py" $?
#
# run_busco ${OUT_FASTA} ${THREADS}
# rm -r busco busco_downloads
#
# busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
#     tee ${OUT_NAME}.csv
#
# # Keep the uncompressed version for output in main directory
# gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
# mv ${OUT_FASTA}.gz ${TARGET}/
#
# cd ..
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz ${TARGET}/
#
# rm -r ${OUT_DIR}
