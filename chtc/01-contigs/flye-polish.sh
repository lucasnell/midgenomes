#!/bin/bash


export ASSEMBLY=contigs_necat.fasta
cp /staging/lnell/assemblies/contigs_necat/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz

# Output names:
export OUT_NAME=contigs_necat_flye-polish
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta



. /app/.bashrc
source /staging/lnell/helpers.sh

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Where to send files:
export TARGET=/staging/lnell/assemblies

mkdir ${OUT_DIR}
cd ${OUT_DIR}
mv ../${ASSEMBLY} ./

export LONGREADS=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?


conda activate flye-env

flye --nano-hq ${LONGREADS} \
    --out-dir flye_out \
    --threads ${THREADS} \
    --iterations 3 \
    --polish-target ${ASSEMBLY}
check_exit_status "flye" $?
conda deactivate

rm ${LONGREADS}

cp flye_out/polished_3.fasta ${OUT_FASTA}
check_exit_status "cp fasta" $?

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads busco

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

