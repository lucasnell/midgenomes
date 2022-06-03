#!/bin/bash

#'
#' Combine assemblies using quickmerge.
#' Requires argument for which assemblies to combine.
#' Also requires argument for the name of the output.
#'
#' Note that the merged assembly size and completeness will be more like the
#' first assembly provided.
#'


# export ASSEMBLY1=contigs_next_ont-polish_nextpolish.fasta
# export ASSEMBLY2=contigs_smart_ont-polish_nextpolish_purgedups.fasta
# export OUT_NAME=contigs_merged_next_smart.fasta
# export SAVE_OUT=1


export ASSEMBLY1=$1
export ASSEMBLY2=$2
export OUT_NAME=$3
# 2 for saving all output, 1 for saving just FASTA, 0 for no saving
export SAVE_OUT=$4


if [[ "${ASSEMBLY1}" != *.fasta ]]; then
    echo "ERROR: Input assembly files must end in *.fasta." 1>&2
    exit 1
fi
if [[ "${ASSEMBLY2}" != *.fasta ]]; then
    echo "ERROR: Input assembly files must end in *.fasta." 1>&2
    exit 1
fi
if [[ "${OUT_NAME}" != *.fasta ]]; then
    echo "ERROR: Output assembly file must end in *.fasta." 1>&2
    exit 1
fi
if ! [[ "${SAVE_OUT}" == "0" ]] && ! [[ "${SAVE_OUT}" == "1" ]] && ! [[ "${SAVE_OUT}" == "2" ]]; then
   echo "ERROR: The 6th input should be 0, 1, or 2. Yours is '${SAVE_OUT}'." 1>&2
   exit 1
fi



. /app/.bashrc
conda activate assembly-env
source /staging/lnell/helpers.sh

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_DIR=${OUT_NAME/.fasta/}
export OUT_FASTA=${OUT_NAME}

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY1}.gz ./ && gunzip ${ASSEMBLY1}.gz
check_exit_status "cp genome 1" $?

cp /staging/lnell/assemblies/${ASSEMBLY2}.gz ./ && gunzip ${ASSEMBLY2}.gz
check_exit_status "cp genome 2" $?


# Make both assemblies not have spaces in sequence names.
sed -i '/^>/{s/ .*//}' ${ASSEMBLY1}
sed -i '/^>/{s/ .*//}' ${ASSEMBLY2}



# Set some parameters for quickmerge:
# Default is 0, but they recommend near the N50 of the first assembly.
export L_CUTOFF=$(summ-scaffs.py ${ASSEMBLY1} | grep "N50" | sed 's/.* //')
# Default is 5000, but they recommend higher
export ML_CUTOFF=10000


merge_wrapper.py -pre ${OUT_DIR} \
    -l ${L_CUTOFF} \
    -ml ${ML_CUTOFF} \
    ${ASSEMBLY1} ${ASSEMBLY2} \
    > >(tee -a quickmerge.out) \
    2> >(tee -a quickmerge.err >&2)

mv merged_${OUT_DIR}.fasta ${OUT_FASTA}


# usage: merge_wrapper.py [-h] [-pre PREFIX] [-hco HCO] [-c C] [-l LENGTH_CUTOFF]
#                              [--no_nucmer] [--no_delta] [--stop_after_nucmer]
#                              [--stop_after_df] [-ml MERGING_LENGTH_CUTOFF]
#                              [--clean_only] hybrid_assembly_fasta self_assembly_fasta
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



# ===================*
# Summarize and move output
# ===================*


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

rm ${ASSEMBLY1} ${ASSEMBLY2}

cd ..

if (( SAVE_OUT == 2 )); then
    tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ${OUT_DIR}.tar.gz ${TARGET}/
fi

rm -r ${OUT_DIR}

