#!/bin/bash

. /app/.bashrc
conda activate assembly-env
source /staging/lnell/helpers.sh

export THREADS=32

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=contigs_necat
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

export LONGREADS=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?

echo $(pwd)/${LONGREADS} > reads.txt


# Make configuration file
#   (Could also run this: `$ necat config tany_config.txt`, then adjust.)
cat <<EOF > tany_config.txt
PROJECT=tany
ONT_READ_LIST=reads.txt
GENOME_SIZE=100000000
THREADS=${THREADS}
MIN_READ_LENGTH=3000
PREP_OUTPUT_COVERAGE=100
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=40
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
EOF


# Correct reads:
necat correct tany_config.txt
check_exit_status "necat correct" $?

# Corrected reads are in the files
# ./tany/1-consensus/cns_iter${NUM_ITER}/cns.fasta

# The longest 40x (CNS_OUTPUT_COVERAGE) corrected reads are extracted for
# assembly, which are in the file
# ./tany/1-consensus/cns_final.fasta

# Assemble reads:
necat assemble tany_config.txt
check_exit_status "necat assemble" $?

# The assembled contigs are in the file
# ./tany/4-fsa/contigs.fasta.

# Bridge contigs:
necat bridge tany_config.txt
check_exit_status "necat bridge" $?

# The bridged contigs are in the file
# ./tany/6-bridge_contigs/bridged_contigs.fasta.

# If POLISH_CONTIGS is set, the pipeline uses the corrected reads to polish
# the bridged contigs. The polished contigs are in the file
# ./tany/6-bridge_contigs/polished_contigs.fasta

cp ./tany/6-bridge_contigs/polished_contigs.fasta ${OUT_FASTA}

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
