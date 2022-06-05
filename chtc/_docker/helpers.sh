#!/bin/bash

#'
#' This script has some two main helper functions used many times throughout.
#'


#' Check previous command's exit status.
#' If != 0, then archive current dir and exit.
#' If the object `TARGET` is defined, it sends the tar file there.
#' The tar file is sent to `/staging/lnell/` otherwise
#' Usage:
#' OperationX ... options ...
#' check_exit_status "OperationX" $?
#'
check_exit_status () {
    if [ "$2" != "0" ]; then
        echo "Step $1 failed with exit status $2" 1>&2
        # Only do exiting if inside a non-interactive job:
        if [ -z "$PS1" ]; then
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
        else
            break 2> /dev/null
            return 0
        fi
    fi
    if [[ "$1" != "null" ]]; then
        echo "Checked step $1"
    fi
    return 0
}



#' Run `busco` on a FASTA file.
#' NOTE: Must be run inside an system containing a `busco-env` conda
#' environment that contains `busco`.
#' Produces directories `busco` and `busco_downloads`, plus the file `busco.out`
#' Usage:
#' run_busco ${OUT_FASTA} ${THREADS}
run_busco () {
    conda activate busco-env
    busco \
        -m genome \
        -l diptera_odb10 \
        -i $1 \
        -o busco \
        --cpu $2 | \
        tee busco.out
    check_exit_status "busco" $?
    conda deactivate
    return 0
}

