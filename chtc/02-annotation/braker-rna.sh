#!/bin/bash

#'
#' Use BRAKER2 for ".. prediction of protein coding gene structures"
#' using RNAseq alignments (pipeline B)
#'
#' Requires arguments for directory containing assembly, output naming,
#' RNA read .tar file(s) or BAM alignment file(s), and assembly name.
#' Optionally takes directory containing read file(s).
#'
#' Usage:
#' braker-rna.sh [-i INPUT_READS_DIR -o OUTPUT_LOC] -p OUTPUT_PREFIX \
#'     -r READS_TAR_FILE_1 -r READS_TAR_FILE_2 ... -r READS_TAR_FILE_N \
#'     ASSEMBLY
#'
#' Options:
#'   -i Directory containing the input read tar files.
#'      Defaults to `/staging/lnell/ill/rna`.
#'   -o Folder for all output. Defaults to `/staging/lnell/annotations`.
#'   -p Prefix for output directory name.
#'      Final output will be `${OUT_PREFIX}_braker_rna.tar.gz`.
#'   -r tar file containing Illumina reads to use, or BAM file containing
#'      alignments from Illumina RNAseq reads.
#'      Must be inside folder given by `-i` arg and end with `.tar` or `.bam`.
#'      For read tar files, this script assumes that only reads are inside the
#'      tar file and that when the file names are sorted alphabetically,
#'      the first file is #1 of pair, and the second is #2 of pair.
#'      For BAM files, it's assumed they're already sorted.
#'      If you want to specify multiple read-tar / BAM files, pass this argument
#'      multiple times.
#'
#'
#' For Tanytarsus gracilentus assembly, the command was...
#'
#' braker-rna.sh -p Tgraci \
#'      -r trimmed_TanyAdult_S1.tar -r trimmed_TanyJuven_S2.tar \
#'      /staging/lnell/annotations/Tgraci_contigs_masked.fasta.gz
#'
#'
#' For Parochlus steinenii assembly, the command was...
#'
#' braker-rna.sh -p Pstein -r trimmed_Pstein_RNA_SRR3951283.tar \
#'      -r trimmed_Pstein_RNA_SRR3951284.tar \
#'      -r trimmed_Pstein_RNA_SRR3951285.tar \
#'      /staging/lnell/annotations/Pstein_contigs_masked.fasta.gz
#'
#'
#' For Culicoides sonorensis assembly, the command was...
#'
#' braker-rna.sh -p Csonor -i /staging/lnell/Cson-RNA \
#'      -r trimmed_ERR637904.tar \
#'      -r trimmed_ERR637905.tar \
#'      -r trimmed_ERR637906.tar \
#'      -r trimmed_ERR637907.tar \
#'      -r trimmed_ERR637908.tar \
#'      -r trimmed_ERR637909.tar \
#'      -r trimmed_ERR637910.tar \
#'      -r trimmed_ERR637911.tar \
#'      -r trimmed_ERR637912.tar \
#'      -r trimmed_ERR637913.tar \
#'      -r trimmed_ERR637914.tar \
#'      /staging/lnell/annotations/Csonor_contigs_masked.fasta.gz
#'



#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================



export READS_LOC=/staging/lnell/ill/rna
export OUTPUT_LOC=/staging/lnell/annotations

unset -v OUT_PREFIX
unset -v TARS_BAMS


while getopts ":i:o:p:r:" opt; do
    case $opt in
        i)
            READS_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        o)
            OUTPUT_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        p)
            OUT_PREFIX="$OPTARG"
            if [[ "$OUT_PREFIX" == */* ]]; then
                echo "ERROR: -p arg cannot contain '/'." 1>&2
                exit 1
            fi
            ;;
        r)
            TARS_BAMS+=("$OPTARG")
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

if [ ! -d "${READS_LOC}" ]; then
    echo "ERROR: Reads directory ('${READS_LOC}') does not exist." 1>&2
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


if [ -z "$OUT_PREFIX" ]; then echo "Missing -p argument." 1>&2; exit 1; fi
if [ -z "$TARS_BAMS" ]; then echo "Missing -r argument(s)." 1>&2; exit 1; fi

export OUT_PREFIX
export TARS_BAMS


tb_status=0
for tb_file in ${TARS_BAMS[@]}; do
    if [ ! -f "${READS_LOC}/${tb_file}" ]; then
        echo "ERROR: '${READS_LOC}/${tb_file}' does not exist." 1>&2
        tb_status=1
    fi
    if ! [[ "$tb_file" =~ (.tar|.bam)$ ]]; then
        echo -n "ERROR: All inputs to -r must end in *.tar or *.bam. " 1>&2
        echo "One of yours is '${tb_file}'." 1>&2
        tb_status=1
    fi
done
if (( tb_status == 1 )); then exit 1; fi
unset tb_status




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

cp ${ASSEMBLY_FULL_PATH} ./
check_exit_status "cp genome" $?
export ASSEMBLY=$(basename "${ASSEMBLY_FULL_PATH}")
if [[ "${ASSEMBLY}" == *.gz ]]; then
    gunzip ${ASSEMBLY}
    check_exit_status "gunzip genome" $?
    ASSEMBLY=${ASSEMBLY%.gz}
fi

export OUT_DIR=${OUT_PREFIX}_braker_rna



#' ===========================================================================
#' ===========================================================================
#'
#' Align RNAseq reads to assembly
#'
#' ===========================================================================
#' ===========================================================================

conda activate main-env


hisat2-build -q ${ASSEMBLY} ${OUT_DIR}_idx
check_exit_status "hisat2-build" $?

export BAM_ARRAY

for tb_file in ${TARS_BAMS[@]}; do
    if [[ "$tb_file" == *.bam ]]; then
        cp ${READS_LOC}/${tb_file} ./
        check_exit_status "cp bam file" $?
        BAM_ARRAY+=("$tb_file")
        continue
    fi
    BAM_FILE=${tb_file%.tar}.bam
    #' Extract individual read names from tar file:
    READS1=$(read_tar_name ${READS_LOC}/${tb_file} 1)
    READS2=$(read_tar_name ${READS_LOC}/${tb_file} 2)
    # Move reads here:
    tar -xf ${READS_LOC}/${tb_file} -C ./
    check_exit_status "extract reads tar file" $?
    if [[ "${READS1}" == *.gz ]]; then
        gunzip ${READS1}
        check_exit_status "gunzip $READS1" $?
        READS1=${READS1%.gz}
    fi
    if [[ "${READS2}" == *.gz ]]; then
        gunzip ${READS2}
        check_exit_status "gunzip $READS2" $?
        READS2=${READS2%.gz}
    fi
    # (Pipe to samtools is to convert it to an aligned BAM file)
    hisat2 -x ${OUT_DIR}_idx \
        -1 ${READS1} -2 ${READS2} \
        -k 3 -p ${THREADS} \
        | samtools sort -O bam - \
        > ${BAM_FILE}
    check_exit_status "hisat2 - $tb_file" $?
    cp ${BAM_FILE} ${OUTPUT_LOC}/
    rm ${READS1} ${READS2}
    BAM_ARRAY+=("$BAM_FILE")
    unset -v BAM_FILE READS1 READS2
done



rm ${OUT_DIR}_idx*


conda deactivate


export BAM_LIST=$(IFS=','; echo "${BAM_ARRAY[*]}")




#' ===========================================================================
#' ===========================================================================
#'
#' Run BRAKER using RNAseq BAM files
#'
#' ===========================================================================
#' ===========================================================================


conda activate annotate-env


braker.pl --genome=${ASSEMBLY} --cores=${THREADS} \
    --softmasking --bam=${BAM_LIST}

mv braker ${OUT_DIR}

# Saving output:
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${OUTPUT_LOC}/

cd ..
rm -r working

