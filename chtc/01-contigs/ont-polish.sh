#!/bin/bash

#'
#' Polish assembly using ONT reads with Racon and medaka.
#' Requires argument for which assembly to polish.
#'

export ASSEMBLY=$1

if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi


. /app/.bashrc
conda activate assembly-env
source /staging/lnell/helpers.sh

export THREADS=16

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_ont-polish
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

export LONGREADS=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?


# ===================*
# Racon
# ===================*


export ALIGN=mm2_align.paf
export NEW_ASSEMBLY=${ASSEMBLY}

# 5 rounds of polishing using Racon:
export npr=5
for i in $(seq 1 ${npr}); do
    minimap2 -x map-ont -t $((THREADS - 1)) -K 1G \
        ${NEW_ASSEMBLY} ${LONGREADS} \
        > ${ALIGN}
    # These parameters below are for compatibility with medaka:
    # -m 8 -x -6 -g -8 -w 500
    racon -m 8 -x -6 -g -8 -w 500 -t ${THREADS} \
        ${LONGREADS} ${ALIGN} ${NEW_ASSEMBLY} \
        > racon_${i}.fasta
    NEW_ASSEMBLY=racon_${i}.fasta
done

rm ${ALIGN}


# Summarize consensus after last round:

summ-scaffs.py racon_${npr}.fasta | tee racon_contigs_summary.out
check_exit_status "summ-scaffs.py (racon)" $?

conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i racon_${npr}.fasta \
    -o busco \
    --cpu ${THREADS} | \
    tee racon_busco.out
check_exit_status "busco (racon)" $?
conda deactivate

rm -r busco_downloads busco

busco_seq_summary_csv racon_contigs_summary.out racon_busco.out \
    ${OUT_NAME}_racon | \
    tee ${OUT_NAME}_racon.csv



# ===================*
# Medaka
# ===================*

#' For the model to use (`-m` argument), the "r941_min_" should remain
#' constant even if you use a later version of guppy because this specifies
#' the pore and device.
medaka_consensus -i ${LONGREADS} -d racon_${npr}.fasta -o medaka -t ${THREADS} \
    -m r941_min_hac_g507

cp ./medaka/consensus.fasta ${OUT_FASTA}

# Summarize medaka consensus:

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS} | \
    tee busco.out
check_exit_status "busco" $?
conda deactivate

rm -r busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv


# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

rm ${LONGREADS}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

