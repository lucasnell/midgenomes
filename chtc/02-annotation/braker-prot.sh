#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using proteins from a database (pipeline C)
#'
#' Requires arguments for directory containing assembly, output naming,
#' and assembly name.
#'
#' Usage:
#' braker-prot.sh -a ASSEMBLY_DIR -o OUTPUT_NAME ASSEMBLY
#'
#' Options:
#'   -a Directory containing the input assembly file.
#'   -o Output directory name.
#'
#'
#' For Tanytarsus gracilentus assembly, the command was...
#'
#' braker-prot.sh -a /staging/lnell/annotation -o tany_braker_prot \
#'      tany_contigs_masked.fasta.gz
#'
#'
#' For Parochlus steinenii assembly, the command was...
#'
#' braker-prot.sh -a /staging/lnell/annotation -o Pstein_braker_prot \
#'      Pstein_contigs_masked.fasta.gz
#'






#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================



unset -v ASSEMBLY_LOC
unset -v OUT_DIR


while getopts ":a:o:" opt; do
    case $opt in
        a)
            ASSEMBLY_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            if [ ! -f "${ASSEMBLY_LOC}" ]; then
                echo "ERROR: '${ASSEMBLY_LOC}' does not exist." 1>&2
                exit 1
            fi
            ;;
        o)
            OUT_DIR="$OPTARG"
            if [[ "$OUT_DIR" == */* ]]; then
                echo "ERROR: -o arg cannot contain '/'." 1>&2
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
if ! [[ "$ASSEMBLY" =~ (.fasta|.fa|.fasta.gz|.fa.gz)$ ]]; then
    echo -n "ERROR: Assembly must end in *.fasta, *.fa, *.fasta.gz, or *.fa.gz. " 1>&2
    echo "Yours is '${ASSEMBLY}'." 1>&2
    exit 1
fi


if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi


if [ -z "$ASSEMBLY_LOC" ]; then echo "Missing -a argument." 1>&2; exit 1; fi
if [ -z "$OUT_DIR" ]; then echo "Missing -o argument." 1>&2; exit 1; fi

export ASSEMBLY_LOC
export OUT_DIR

if [ ! -f ${ASSEMBLY_LOC}/${ASSEMBLY} ]; then
    echo "ERROR: '${ASSEMBLY_LOC}/${ASSEMBLY}' does not exist." 1>&2
    exit 1
fi



#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation
if [ ! -f ${TARGET} ]; then
    echo "INTERNAL ERROR: '${TARGET}' (TARGET object) does not exist." 1>&2
    exit 1
fi

mkdir working
cd working


cp ${ASSEMBLY_LOC}/${ASSEMBLY} ./
check_exit_status "cp genome" $?
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
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working

