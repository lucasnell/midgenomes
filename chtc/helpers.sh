#!/bin/bash

#'
#' This script has some helper functions used throughout.
#'


#' Check previous command's exit status.
#' If != 0, then archive current dir and exit.
#' If the object `TARGET` is defined, it sends the tar file there.
#' The tar file is sent to `/staging/lnell/` otherwise
#' Usage:
#' check_exit_status "NameOfOperation" $?
check_exit_status () {
    if [ ! "$2" -eq "0" ]; then
        echo "Step $1 failed with exit status $2" 1>&2
        local wd=${PWD##*/}
        cd ..
        tar -czf ERROR_${wd}.tar.gz ${wd}
        local SEND_LOC=/staging/lnell/
        if [ -n "$TARGET" ]; then
            SEND_LOC=${TARGET}
        fi
        mv ERROR_${wd}.tar.gz ${SEND_LOC}
        rm -r ${wd}
        exit $2
    fi
    echo "Checked step $1"
}

#' Check on BAM file with `bamtools stats`, check status of this call,
#' then output the file name.
#' Usage:
#' call_bam_stats ${BAM_FILE} "description"
call_bam_stats () {
    local B=$1
    local S=${B/.bam/.stats}
    bamtools stats -in $B | tee $S
    check_exit_status "bamtools stats $2" $?
    echo -e "FILE:" $B "\n**********************************************"
}

#'
#' Organize output from `summ-scaffs.py` and `busco` into a string
#' that can be used in a CSV file.
#' Usage:
#' busco_seq_summary_csv ${SUMM_SCAFFS_OUT} ${BUSCO_OUT} ${SEQS_DESCRIPTOR}
busco_seq_summary_csv () {
    local OUT=""
    # Header:
    OUT+="file,size,scaffs,N50,min,max,total_N,"
    OUT+="BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n\n"
    # Output file identifier:
    OUT+="${3},"
    # Output from summ-scaffs.py:
    OUT+="$(grep "size" ${1} | sed 's/.* //'),"
    OUT+="$(grep "scaffolds$" ${1} | sed -r 's/\ .+//'),"
    OUT+="$(grep "N50" ${1} | sed 's/.* //'),"
    OUT+="$(grep "min" ${1} | sed 's/.* //'),"
    OUT+="$(grep "max" ${1} | sed 's/.* //'),"
    OUT+="$(grep "total\ N" ${1} | sed 's/.* //'),"
    # Output from BUSCO:
    OUT+="$(grep '|' ${2} | grep '(C)' | tr -d '|' | xargs | sed 's/ .*//'),"
    OUT+="$(grep '|' ${2} | grep '(S)' | tr -d '|' | xargs | sed 's/ .*//'),"
    OUT+="$(grep '|' ${2} | grep '(D)' | tr -d '|' | xargs | sed 's/ .*//'),"
    OUT+="$(grep '|' ${2} | grep '(F)' | tr -d '|' | xargs | sed 's/ .*//'),"
    OUT+="$(grep '|' ${2} | grep '(M)' | tr -d '|' | xargs | sed 's/ .*//'),"
    OUT+="$(grep '|' ${2} | grep 'Total BUSCO' | tr -d '|' | xargs | sed 's/ .*//')\n"
    # Use printf bc it provides actual newlines
    printf ${OUT}
}



