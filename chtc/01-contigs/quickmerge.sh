#!/bin/bash

#'
#' Combine assemblies using quickmerge.
#' Requires argument for which assemblies to combine.
#' Also requires argument for the name of the output.
#'
#' Note that the merged assembly size and completeness will be more like the
#' first assembly provided.
#'

export ASSEMBLY1=$1
export ASSEMBLY2=$2
export OUT_NAME=$3

if [[ "${ASSEMBLY1}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi
if [[ "${ASSEMBLY2}" != *.fasta ]]; then
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
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY1}.gz ./ && gunzip ${ASSEMBLY1}.gz
check_exit_status "cp genome 1" $?

cp /staging/lnell/assemblies/${ASSEMBLY2}.gz ./ && gunzip ${ASSEMBLY2}.gz
check_exit_status "cp genome 2" $?


# Also, make both assemblies not have spaces in sequence names.


# try the command 'merge_wrapper.py -h' for detail on options available with this wrapper.
# merge_wrapper.py hybrid_assembly.fasta self_assembly.fasta

# > merge_wrapper.py -h
# usage: merge_wrapper.py [-h] [-pre PREFIX] [-hco HCO] [-c C] [-l LENGTH_CUTOFF] [--no_nucmer] [--no_delta] [--stop_after_nucmer] [--stop_after_df] [-ml MERGING_LENGTH_CUTOFF] [--clean_only]
#                         hybrid_assembly_fasta self_assembly_fasta
#
# run mummer and the merge program.
#
# positional arguments:
#   hybrid_assembly_fasta
#                         the output of a hybrid assembly program such as DBG2OLC
#   self_assembly_fasta   the output of a self assembly program such as PBcR
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -pre PREFIX, --prefix PREFIX
#                         the prefix for all output files
#   -hco HCO, --hco HCO   the quickmerge hco parameter (default=5.0)
#   -c C, --c C           the quickmerge c parameter (default=1.5)
#   -l LENGTH_CUTOFF, --length_cutoff LENGTH_CUTOFF
#                         minimum seed contig length to be merged (default=0)
#   --no_nucmer           skip the nucmer step
#   --no_delta            skip the nucmer and delta-filter steps
#   --stop_after_nucmer   do not perform the delta-filter and merger steps
#   --stop_after_df       do not perform the merger step
#   -ml MERGING_LENGTH_CUTOFF, --merging_length_cutoff MERGING_LENGTH_CUTOFF
#                         set the merging length cutoff necessary for use in quickmerge (default 5000)
#   --clean_only          generate safe FASTA files for merging, but do not merge
#





# ===================*
# Summarize and move output
# ===================*



cp XXXXXXXXXXXXXXXXXXX ${OUT_FASTA}

# Summarize new assembly

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv


# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

rm ${ASSEMBLY1} ${ASSEMBLY2}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

