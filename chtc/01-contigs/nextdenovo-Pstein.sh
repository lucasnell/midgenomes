#!/bin/bash

. /app/.bashrc


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Memory available in GB - should be > 48, ideally >= 64
export MEMORY=$(grep "^Memory = " $_CONDOR_MACHINE_AD | sed 's/Memory\ =\ //')
MEMORY=$(( MEMORY / 1000 ))

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=Pstein_contigs_next
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

export LONGREADS=Pstein_20180703DL1_LabG03.fastq
cp /staging/lnell/ont/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?


echo ${LONGREADS} > input.fofn

# Setting some of the threading and memory options based on
# https://nextdenovo.readthedocs.io/en/latest/FAQ.html#how-to-optimize-parallel-computing-parameters

TOTAL_INPUT_BASES=20
# Below, it can be 32-64, so I chose 40
PJ=$(python -c "print(round(${MEMORY} / 40))")
PC=$(python -c "print(round(${MEMORY} / (${TOTAL_INPUT_BASES} * 1.2/4)))")
SOm=$(python -c "print(round(${TOTAL_INPUT_BASES} * 1.2/4))")
SOt=$(python -c "print(round(${THREADS} / ${PC}))")
CO=${SOt}
MOR=$(python -c "print(round(${THREADS} / ${PJ}))")
MOC=${MOR}


cat << EOF > pstein.cfg
[General]
job_type = local
job_prefix = ${OUT_NAME}
task = all
rewrite = yes
deltmp = yes
parallel_jobs = ${PJ}
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = workdir

[correct_option]
read_cutoff = 1k
genome_size = 144M
sort_options = -m ${SOm}G -t ${SOt}
minimap2_options_raw = -t ${MOR}
pa_correction = ${PC}
correction_options = -p ${CO}

[assemble_option]
minimap2_options_cns = -t ${MOC}
nextgraph_options = -a 1
EOF



nextDenovo pstein.cfg
check_exit_status "nextDenovo" $?


#' Output:
#'
#' workdir/03.ctg_graph/nd.asm.fasta
#' - Contigs with fasta format, the fasta header includes ID, type, length,
#'   node count, a consecutive lowercase region in the sequence implies a
#'   weak connection, and a low quality base is marked with a single
#'   lowercase base.
#' workdir/03.ctg_graph/nd.asm.fasta.stat
#' - Some basic statistical information (N10-N90, Total size et al.).
#'

cp workdir/03.ctg_graph/nd.asm.fasta ${OUT_FASTA}
check_exit_status "rename FASTA" $?

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads

pretty-csv.py -s contigs_summary.out -b busco.out ${OUT_NAME} \
    | tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

rm ${LONGREADS}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
