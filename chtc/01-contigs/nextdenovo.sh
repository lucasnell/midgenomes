#!/bin/bash

. /app/.bashrc
source /staging/lnell/helpers.sh

export THREADS=32

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=contigs_next
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

export LONGREADS=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
check_exit_status "cp guppy" $?


echo ${LONGREADS} > input.fofn

# Could also the code below, then adjust manually.
# cp /opt/NextDenovo/doc/run.cfg ./tany.cfg

cat <<EOF > tany.cfg
[General]
job_type = local
job_prefix = ${OUT_NAME}
task = all
rewrite = yes
deltmp = yes
parallel_jobs = ${THREADS}
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = workdir

[correct_option]
read_cutoff = 1k
genome_size = 100M
sort_options = -m 20g -t 15
minimap2_options_raw = -t 8
pa_correction = 6 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -t 8
nextgraph_options = -a 1

# see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters
EOF


nextDenovo run.cfg
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
