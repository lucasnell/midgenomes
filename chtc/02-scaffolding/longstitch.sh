#!/bin/bash


export THREADS=16

# Input FASTA is given by submit file (NO EXTENSION SHOULD BE GIVEN):
export GENOME=$1
# Values for arguments k_ntLink and w, respectively.
export K=$2
export W=$3
# 1 for saving output, 0 for not
export SAVE_OUT=$4

export TARGET=/staging/lnell/assemblies

export init_code=0
if [ ! -f ${TARGET}/${GENOME}.fasta.gz ]; then
    echo -e "\n\nERROR: ${TARGET}/${GENOME}.fasta.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    init_code=1
fi
if ! [[ $K =~ ^[0-9]+$ ]]; then
    echo "ERROR: Second arg is not an integer! " 1>&2
    init_code=1
fi
if ! [[ $W =~ ^[0-9]+$ ]]; then
    echo "ERROR: Third arg is not an integer! " 1>&2
    init_code=1
fi
if ! [[ "${SAVE_OUT}" == "0" ]] && ! [[ "${SAVE_OUT}" == "1" ]]; then
   echo "ERROR: The 6th input should be 0 or 1. Yours is '${SAVE_OUT}'." 1>&2
   init_code=1
fi
if (( init_code != 0 )); then
    exit 1
fi



. /app/.bashrc
source /staging/lnell/helpers.sh
conda activate main-env




# The input file dictates the output name:
export OUT_SUFFIX=_longstitch

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME}${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta




mkdir ${OUT_DIR}

cd ${OUT_DIR}

cp ${TARGET}/${GENOME}.fasta.gz ./ && gunzip ${GENOME}.fasta.gz
# LongStitch requires *.fa ending
mv ${GENOME}.fasta ${GENOME}.fa

export FASTQ=basecalls_guppy-5.0.11
cp /staging/lnell/${FASTQ}.fastq.gz ./
# LongStitch requires *.fq.gz ending
mv ${FASTQ}.fastq.gz ${FASTQ}.fq.gz


# longstitch run
conda activate longstitch-env
longstitch tigmint-ntLink-arks \
    draft=${GENOME} reads=${FASTQ} t=${THREADS} G=100000000 \
    k_ntLink=${K} w=${W}
check_exit_status "longstitch" $?
conda deactivate

rm -rf ${GENOME}.fa ${FASTQ}.fq.gz

if [ ! -f *-arks.longstitch-scaffolds.fa ]; then
    echo -e "\n\nERROR: LongStitch failed to produce output." 1>&2
    echo -e "Exiting...\n" 1>&2
    cd ..
    rm -r ./${OUT_DIR}
    exit 1
fi



# Moving the final FASTA to staging
# Gets filename of final scaffolds FASTA file
export LS_FASTA=$(basename -- $(readlink -f *-arks.longstitch-scaffolds.fa))
cp $LS_FASTA ${OUT_FASTA}

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads busco

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_DIR} \
    | tee ${OUT_DIR}.csv


if (( SAVE_OUT == 0 )); then
    cd ..
    rm -r ${OUT_DIR}
    exit 0
fi


# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


