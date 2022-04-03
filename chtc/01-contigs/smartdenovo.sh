#!/bin/bash

. /app/.bashrc
source /staging/lnell/helpers.sh

export THREADS=32

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=contigs_smart
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/${LONGREADS} ./
check_exit_status "cp guppy" $?

conda activate main-env
seqtk seq -A -l 80 ${LONGREADS} > reads.fa
check_exit_status "seqtk reads" $?
rm ${LONGREADS}

conda activate assembly-env

smartdenovo.pl -p tany -c 1 -t ${THREADS} reads.fa > tany.mak
make -f tany.mak

# the raw unitigs are reported in file
# tany.lay.utg

# consensus unitigs in file
# tany.cns


cp tany.dmo.cns ${OUT_FASTA}

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

rm reads.fa tany.fa.gz
# Normally I'd keep some intermediates, but these total > 80 GB
rm tany.dmo*

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

