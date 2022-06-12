#!/bin/bash

#'
#' This script has some two main helper functions used many times throughout.
#'


#' Check previous command's exit status.
#' If != 0, then...
#' ... non-interactive session: send error message, print any *.log or
#'                              *.err files, and exit.
#' ... interactive session: just break and ignore errors if the context
#'                          isn't sensible.
#'
#' Usage:
#' OperationX ... options ...
#' check_exit_status "OperationX" $?
#'
check_exit_status () {
    if [ "$2" != "0" ]; then
        echo "Step $1 failed with exit status $2" 1>&2
        local LOGS=$(ls *.err *.out 2> /dev/null)
        for f in $LOGS; do
            echo -e "\n\nFILE: " ${f} "\n\n" 1>&2
            cat ${f} 1>&2
        done
        out_loc=$(readlink -f /proc/self/fd/0)
        if [ "$out_loc" != "/dev/null" ]; then
            # (interactive shell)
            break 2> /dev/null
        else
            # (non-interactive shell)
            cd $_CONDOR_SCRATCH_DIR
            # Remove everything except for the files / folders you start with:
            TO_RM=$(ls -I tmp -I var -I docker_stderror)
            rm ${TO_RM} 2> /dev/null
            exit $2
        fi
        return 0
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

