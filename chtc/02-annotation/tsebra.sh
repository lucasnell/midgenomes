#!/bin/bash


#'
#' Combine 2 BRAKER annotations (RNAseq + OrthoDB proteins) using TSEBRA,
#' and run BUSCO on resulting CDS.
#' The RNAseq annotation will be preferred because we're assuming the
#' proteins from OrthoDB are generally not from closely related species.
#'
#' Requires arguments for both annotations, output naming, and assembly.
#' It also takes an optional argument for the directory where output should go.
#'
#' Usage:
#' tsebra.sh [-o OUTPUT_LOC] -r RNASEQ_ANNOTATION -t ORTHO_ANNOTATION \
#'    -p OUT_PREFIX ASSEMBLY
#'
#' Options:
#'   -o Folder for all output. Defaults to `/staging/lnell/annotations`.
#'   -r Full path to compressed folder containing BRAKER annotation using
#'      RNAseq data. Must end in `.tar.gz`.
#'      It's assumed that the folder name is the same as the compressed tarball
#'      minus the `.tar.gz`.
#'      Inside the folder, there should be files `augustus.hints.gtf`
#'      and `hintsfile.gff`.
#'   -t Compressed folder containing BRAKER annotation using OrthoDB proteins.
#'      The notes for argument `-r` also apply here.
#'   -p Prefix for all output. Final outputs will be
#'      `${OUT_PREFIX}_tsebra.gff3.gz`,
#'      `${OUT_PREFIX}_cds.fasta.gz`,
#'      `${OUT_PREFIX}_proteins.faa.gz`, and
#'      `${OUT_PREFIX}_cds_busco.tar.gz`.
#'
#'



#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotations

unset -v RNASEQ_ANNOT_FULL_PATH
unset -v ORTHO_ANNOT_FULL_PATH
unset -v OUT_PREFIX


while getopts ":o:r:t:p:" opt; do
    case $opt in
        o)
            OUTPUT_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        r)
            RNASEQ_ANNOT_FULL_PATH="$OPTARG"
            ;;
        t)
            ORTHO_ANNOT_FULL_PATH="$OPTARG"
            ;;
        p)
            OUT_PREFIX="$OPTARG"
            if [[ "$OUT_PREFIX" == */* ]]; then
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
fi
if ! [[ "${ASSEMBLY_FULL_PATH}" =~ (.fasta|.fa|.fasta.gz|.fa.gz)$ ]]; then
    echo -n "ERROR: Assembly must end in '.fasta', '.fa', '.fasta.gz', or '.fa.gz'. " 1>&2
    echo "Yours is '${ASSEMBLY_FULL_PATH}'." 1>&2
    exit 1
fi

if [ ! -d "${OUTPUT_LOC}" ]; then
    echo "ERROR: Output directory ('${OUTPUT_LOC}') does not exist." 1>&2
    exit 1
fi

if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi


if [ -z "$RNASEQ_ANNOT_FULL_PATH" ]; then echo "Missing -r argument(s)." 1>&2; exit 1; fi
if [ -z "$ORTHO_ANNOT_FULL_PATH" ]; then echo "Missing -t argument(s)." 1>&2; exit 1; fi
if [ -z "$OUT_PREFIX" ]; then echo "Missing -p argument." 1>&2; exit 1; fi


export RNASEQ_ANNOT_FULL_PATH
export ORTHO_ANNOT_FULL_PATH
export OUT_PREFIX

for anno in "${RNASEQ_ANNOT_FULL_PATH}" "${ORTHO_ANNOT_FULL_PATH}"; do
    if [ ! -f "${anno}" ]; then
        echo "ERROR: Your annotation '${anno}' does not exist." 1>&2
        exit 1
    fi
    if [[ "${anno}" != *.tar.gz ]]; then
        echo "ERROR: Annotations must end in '.tar.gz'. Yours is '${anno}'." 1>&2
        exit 1
    fi
done







#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================



export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc



#' Prefix (everything but extension) to eventual gff3 file:
export OUT_GFF_PREFIX=${OUT_PREFIX}_tsebra
#' FASTA file containing coding DNA sequences (CDS):
export OUT_CDS_FASTA=${OUT_PREFIX}_cds.fasta
#' FASTA file (.faa to show it's AA sequences) containing proteins:
export OUT_PROT_FASTA=${OUT_PREFIX}_proteins.faa

#' Output directory for BUSCO run on CDS:
export BUSCO_OUT_DIR=${OUT_PREFIX}_cds_busco


mkdir working
cd working

#' Small function to perform if something fails
disgraceful_exit () {
    cd $_CONDOR_SCRATCH_DIR
    rm -r working
    exit 1
}





#' ===========================================================================
#' ===========================================================================
#'
#' Copy over BRAKER1 and BRAKER2 gene predictions:
#'
#' ===========================================================================
#' ===========================================================================


tar -xzf ${RNASEQ_ANNOT_FULL_PATH} -C ./
check_exit_status "cp, extract rnaseq annotation" $?
export RNASEQ_ANNOT=$(basename ${RNASEQ_ANNOT_FULL_PATH})
RNASEQ_ANNOT=${RNASEQ_ANNOT%.tar.gz}

tar -xzf ${ORTHO_ANNOT_FULL_PATH} -C ./
check_exit_status "cp, extract ortho annotation" $?
export ORTHO_ANNOT=$(basename ${ORTHO_ANNOT_FULL_PATH})
ORTHO_ANNOT=${ORTHO_ANNOT%.tar.gz}

# ---------------
# Check resulting folders / files:
# ---------------
if [ ! -d ${RNASEQ_ANNOT} ]; then
    echo -n "Decompression and un-tarring RNAseq annotation should " 1>&2
    echo "result in the folder '${RNASEQ_ANNOT}'." 1>&2
    disgraceful_exit
fi
if [ ! -d ${ORTHO_ANNOT} ]; then
    echo -n "Decompression and un-tarring OrthoDB annotation should " 1>&2
    echo "result in the folder '${ORTHO_ANNOT}'." 1>&2
    disgraceful_exit
fi
if [ ! -f ${RNASEQ_ANNOT}/augustus.hints.gtf ] || [ ! -f ${RNASEQ_ANNOT}/hintsfile.gff ]; then
    echo -n "Folder '${RNASEQ_ANNOT}' should contain files " 1>&2
    echo "'augustus.hints.gtf' and 'hintsfile.gff'." 1>&2
    disgraceful_exit
fi
if [ ! -f ${ORTHO_ANNOT}/augustus.hints.gtf ] || [ ! -f ${ORTHO_ANNOT}/hintsfile.gff ]; then
    echo -n "Folder '${ORTHO_ANNOT}' should contain files " 1>&2
    echo "'augustus.hints.gtf' and 'hintsfile.gff'." 1>&2
    disgraceful_exit
fi





#' ===========================================================================
#' ===========================================================================
#'
#' Run TSEBRA
#'
#' ===========================================================================
#' ===========================================================================


conda activate main-env

#' Because the proteins used in the OrthoDB annotations aren't from
#' closely related species, we'll prefer the predictions from RNAseq.
#' Hence the use of `pref_braker1.cfg` below.

tsebra.py \
    -g ${RNASEQ_ANNOT}/augustus.hints.gtf,${ORTHO_ANNOT}/augustus.hints.gtf \
    -c /opt/TSEBRA/config/pref_braker1.cfg \
    -e ${RNASEQ_ANNOT}/hintsfile.gff,${ORTHO_ANNOT}/hintsfile.gff \
    -o ${OUT_GFF_PREFIX}.gtf

#' This fixes IDs in output GTF, which allows it to be converted to GFF3
fix_gtf_ids.py --gtf ${OUT_GFF_PREFIX}.gtf --out ${OUT_GFF_PREFIX}_fixed.gtf
mv ${OUT_GFF_PREFIX}_fixed.gtf ${OUT_GFF_PREFIX}.gtf

conda deactivate


#' ------------------------------------------------------------
#' Create final output files:
#' ------------------------------------------------------------

conda activate annotate-env

gffread -E ${OUT_GFF_PREFIX}.gtf > ${OUT_GFF_PREFIX}.gff3

# We'll only use the gff3 version from now on.
rm ${OUT_GFF_PREFIX}.gtf


cp ${ASSEMBLY_FULL_PATH} ./
check_exit_status "cp assembly" $?
export ASSEMBLY=$(basename "${ASSEMBLY_FULL_PATH}")
if [[ "${ASSEMBLY}" == *.gz ]]; then
    gunzip ${ASSEMBLY}
    check_exit_status "gunzip assembly" $?
    ASSEMBLY=${ASSEMBLY%.gz}
fi


# "spliced CDS for each GFF transcript"
gffread -x ${OUT_CDS_FASTA} -g ${ASSEMBLY} ${OUT_GFF_PREFIX}.gff3
# "protein fasta file with the translation of CDS for each record"
gffread -y ${OUT_PROT_FASTA} -g ${ASSEMBLY} ${OUT_GFF_PREFIX}.gff3

conda deactivate







#' ===========================================================================
#' ===========================================================================
#'
#' Run BUSCO on resulting CDS and
#' move over final outputs
#'
#' ===========================================================================
#' ===========================================================================

run_busco ${OUT_CDS_FASTA} ${THREADS}

mv busco ${BUSCO_OUT_DIR}
mv busco.out ./${BUSCO_OUT_DIR}/busco.stdout

tar -czf ${BUSCO_OUT_DIR}.tar.gz ${BUSCO_OUT_DIR}
mv ${BUSCO_OUT_DIR}.tar.gz ${OUTPUT_LOC}/


#' ----------
#' Move over final outputs:
#' ----------


gzip ${OUT_GFF_PREFIX}.gff3
gzip ${OUT_CDS_FASTA}
gzip ${OUT_PROT_FASTA}

mv ${OUT_GFF_PREFIX}.gff3.gz ${OUT_CDS_FASTA}.gz ${OUT_PROT_FASTA}.gz ${OUTPUT_LOC}/

cd $_CONDOR_SCRATCH_DIR
rm -rf working
exit 0

