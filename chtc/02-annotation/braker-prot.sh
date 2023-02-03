#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using proteins from a database (pipeline C)
#'
#' Requires arguments for output naming and assembly file.
#'
#' Usage:
#' braker-prot.sh [-o OUTPUT_LOC] -p OUTPUT_PREFIX ASSEMBLY
#'
#' Options:
#'   -o Folder for all output. Defaults to `/staging/lnell/annotation`.
#'   -p Prefix for output directory name.
#'      Final output will be `${OUT_PREFIX}_braker_prot.tar.gz`.
#'
#'
#' For Tanytarsus gracilentus assembly, the command was...
#'
#' braker-prot.sh -p Tgraci \
#'      /staging/lnell/annotation/Tgraci_contigs_masked.fasta.gz
#'
#'
#' For Parochlus steinenii assembly, the command was...
#'
#' braker-prot.sh -p Pstein \
#'      /staging/lnell/annotation/Pstein_contigs_masked.fasta.gz
#'
#'
#' For Culicoides sonorensis assembly, the command was...
#'
#' braker-prot.sh -p Csonor \
#'      /staging/lnell/annotation/Csonor_contigs_masked.fasta.gz
#'






#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotation

unset -v OUT_PREFIX


while getopts ":o:p:" opt; do
    case $opt in
        o)
            OUTPUT_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        p)
            OUT_PREFIX="$OPTARG"
            if [[ "${OUT_PREFIX}" == */* ]]; then
                echo "ERROR: -p arg cannot contain '/'." 1>&2
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

export ASSEMBLY_FULL_PATH=${@:$OPTIND:1}
if [ ! -f "${ASSEMBLY_FULL_PATH}" ]; then
    echo "ERROR: Your assembly file ('${ASSEMBLY_FULL_PATH}') does not exist." 1>&2
    exit 1
if ! [[ "${ASSEMBLY_FULL_PATH}" =~ (.fasta|.fa|.fasta.gz|.fa.gz)$ ]]; then
    echo -n "ERROR: Assembly must end in '.fasta', '.fa', '.fasta.gz', or '.fa.gz'. " 1>&2
    echo "Yours is '${ASSEMBLY_FULL_PATH}'." 1>&2
    exit 1
fi
fi

if [ ! -d "${OUTPUT_LOC}" ]; then
    echo "ERROR: Output directory ('${OUTPUT_LOC}') does not exist." 1>&2
    exit 1
fi

if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi


if [ -z "$OUT_PREFIX" ]; then echo "Missing -p argument." 1>&2; exit 1; fi

export OUT_PREFIX





#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc


mkdir working
cd working

export OUT_DIR=${OUT_PREFIX}_braker_prot


cp ${ASSEMBLY_FULL_PATH} ./
check_exit_status "cp genome" $?
export ASSEMBLY=$(basename "${ASSEMBLY_FULL_PATH}")
if [[ "${ASSEMBLY}" == *.gz ]]; then
    gunzip ${ASSEMBLY}
    ASSEMBLY=${ASSEMBLY%.gz}
    check_exit_status "gunzip genome" $?
fi





#' ===========================================================================
#' ===========================================================================
#'
#' Run BRAKER2
#'
#' ===========================================================================
#' ===========================================================================



# Using arthropoda proteins from OrthoDB v10:
wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar -xzf odb10_arthropoda_fasta.tar.gz
cat arthropoda/Rawdata/* > odb10_arthropoda.fasta
rm -r arthropoda odb10_arthropoda_fasta.tar.gz


conda activate annotate-env

braker.pl --genome=${ASSEMBLY} --prot_seq=odb10_arthropoda.fasta \
    --softmasking --cores=${THREADS}
check_exit_status "braker" $?

mv braker ${OUT_DIR}

# Saving output:
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${OUTPUT_LOC}/

cd ..
rm -r working

