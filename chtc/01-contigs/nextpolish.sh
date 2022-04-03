#!/bin/bash

#'
#' Polish assembly using Illumina reads with NextPolish.
#' Requires argument for which assembly to polish.
#'

export ASSEMBLY=$1

if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi


# . /app/.bashrc
# conda activate assembly-env
source /staging/lnell/helpers.sh


export THREADS=16

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_nextpolish
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

export READS1=trimmed_MyKS-19-B_S18_L002_R1_001.fastq
export READS2=trimmed_MyKS-19-B_S18_L002_R2_001.fastq
cp /staging/lnell/dna/trimmed/${READS1}.gz ./ && gunzip ${READS1}.gz
check_exit_status "cp reads 1" $?
cp /staging/lnell/dna/trimmed/${READS2}.gz ./ && gunzip ${READS2}.gz
check_exit_status "cp reads 2" $?


echo ${READS1} ${READS2} > sgs.fofn

echo -e "task = best\ngenome = ${ASSEMBLY}\nsgs_fofn = sgs.fofn" > run.cfg

nextPolish run.cfg
check_exit_status "nextPolish" $?

# Finally polished genome:
# Sequence: /path_to_work_directory/genome.nextpolish.fasta
# Statistics: /path_to_work_directory/genome.nextpolish.fasta.stat

cp genome.nextpolish.fasta ${OUT_FASTA}
check_exit_status "cp output" $?

rm ${READS1} ${READS2}

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

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
