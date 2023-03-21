#!/bin/bash


#'
#' Functional annotations using mantis.
#'
#' Requires arguments for output naming, lineage name, and proteins FASTA.
#' It also takes optional arguments for the directory containing assembly
#' and the directory where output should go.
#'
#' Usage:
#' mantis.sh [-o OUTPUT_LOC -m MANTIS_DL] \
#'            -p OUT_PREFIX -n NCBI_LINEAGE PROTEINS
#'
#' Options:
#'   -o Folder for all output. Defaults to `/staging/lnell/annotations/`.
#'   -m Full path to compressed tar file containing mantis databases that
#'      were downloaded previously. Defaults to
#'      `/staging/lnell/annotations/mantis-downloads.tar.gz`.
#'   -p Prefix for all output.
#'      Final output will be `${OUT_PREFIX}_mantis.tar.gz`.
#'   -n NCBI lineage name. Should be entirely numeric.
#'      Look up from https://www.ncbi.nlm.nih.gov/taxonomy
#'
#'
#' For Tanytarsus gracilentus assembly, the command was...
#' mantis.sh -p Tgraci -n 288803 /staging/lnell/annotations/Tgraci_proteins.faa.gz
#'
#' For Parochlus steinenii assembly, the command was...
#' mantis.sh -p Pstein -n 315571 /staging/lnell/annotations/Pstein_proteins.faa.gz
#'
#' For Chironomus riparius assembly, the command was...
#' mantis.sh -p Cripar -n 315576 /staging/lnell/annotations/Cripar_proteins.faa.gz
#'
#' For Culicoides sonorensis assembly, the command was...
#' mantis.sh -p Csonor -n 179676 /staging/lnell/annotations/Csonor_proteins.faa.gz
#'


export OUTPUT_LOC=/staging/lnell/annotations
export MANTIS_DL_FULL_PATH=/staging/lnell/annotations/mantis-downloads.tar.gz


unset -v OUT_PREFIX
unset -v NCBI_LINEAGE

while getopts ":o:m:p:n:" opt; do
    case $opt in
        o)
            OUTPUT_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        m)
            MANTIS_DL_FULL_PATH="$OPTARG"
            if [[ "${MANTIS_DL_FULL_PATH}" != *.tar.gz ]]; then
                echo "ERROR: -m arg must end in '.tar.gz'." 1>&2
                exit 1
            fi
            ;;
        p)
            OUT_PREFIX="$OPTARG"
            if [[ "OUT_PREFIX" == */* ]]; then
                echo "ERROR: -p arg cannot contain '/'." 1>&2
                exit 1
            fi
            ;;
        n)
            NCBI_LINEAGE="$OPTARG"
            if [[ ! "$NCBI_LINEAGE" =~ ^[0-9]+$ ]]; then
                echo "ERROR: -n arg must be entirely numeric." 1>&2
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

export PROTEINS_FULL_PATH=${@:$OPTIND:1}
if [ ! -f "${PROTEINS_FULL_PATH}" ]; then
    echo "ERROR: Your proteins file ('${PROTEINS_FULL_PATH}') does not exist." 1>&2
    exit 1
fi
if ! [[ "${PROTEINS_FULL_PATH}" =~ (.faa|.faa.gz)$ ]]; then
    echo -n "ERROR: Proteins must end in *.faa, or *.faa.gz. " 1>&2
    echo "Yours is '${PROTEINS_FULL_PATH}'." 1>&2
    exit 1
fi

if [ ! -d "${OUTPUT_LOC}" ]; then
    echo "ERROR: Output directory ('${OUTPUT_LOC}') does not exist." 1>&2
    exit 1
fi
if [ ! -f "${MANTIS_DL_FULL_PATH}" ]; then
    echo "ERROR: '${MANTIS_DL_FULL_PATH}' does not exist." 1>&2
    exit 1
fi

if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi


if [ -z "$OUT_PREFIX" ]; then echo "Missing -p argument." 1>&2; exit 1; fi
if [ -z "$NCBI_LINEAGE" ]; then echo "Missing -n argument." 1>&2; exit 1; fi

export OUT_PREFIX
export NCBI_LINEAGE



#' ===========================================================================
#' ===========================================================================
#'
#' Prepare for mantis run
#'
#' ===========================================================================
#' ===========================================================================


. /app/.bashrc
conda activate annotate-env

export THREADS=$(count_threads)

# Memory available:
export MEMORY=$(grep "^Memory = " $_CONDOR_MACHINE_AD | sed 's/Memory\ =\ //')
# In GB, with 5 GB overhead:
MEMORY=$(( MEMORY / 1000 - 5 ))

export OUT_DIR=${OUT_PREFIX}_mantis

mkdir working
cd working



cp ${PROTEINS_FULL_PATH} ./
check_exit_status "cp proteins" $?
export PROTEINS=$(basename "${PROTEINS_FULL_PATH}")
if [[ "${PROTEINS}" == *.gz ]]; then
    gunzip ${PROTEINS}
    check_exit_status "gunzip proteins" $?
    PROTEINS=${PROTEINS%.gz}
fi

tar -xf ${MANTIS_DL_FULL_PATH} -C ./
check_exit_status "cp, extract databases" $?
if [ ! -d "dbs" ] || [ ! -d "refs" ]; then
    echo -n "ERROR: when un-tarred, '${MANTIS_DL_FULL_PATH}' should create " 1>&2
    echo "'dbs' and 'refs' folders." 1>&2
    cd ..
    rm -r working
    exit 1
fi


#' Setup MANTIS.cfg based on default from
#' https://github.com/PedroMTQ/mantis/blob/a59937d23979372c65188cf2358b69485099dc68/config/MANTIS.cfg
#' downloaded on 14 June 2022.
#' The only differences are that comments are removed and
#' default_ref_folder and resources_folder are set.
cat << EOF > MANTIS.cfg
default_ref_folder=$(pwd)/dbs
resources_folder=$(pwd)/refs
nog_ref=dmnd
nog_weight=0.8
pfam_weight=0.9
EOF

#' ===========================================================================
#' ===========================================================================
#'
#' Run mantis and move output
#'
#' ===========================================================================
#' ===========================================================================

#' Run mantis:
mantis run -mc MANTIS.cfg -i ${PROTEINS} -o ${OUT_DIR} -od ${NCBI_LINEAGE} \
    --verbose_kegg_matrix \
    --cores ${THREADS} --memory ${MEMORY} --hmmer_threads 2

rm -rf refs dbs


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${OUTPUT_LOC}/

cd ..
rm -r working






