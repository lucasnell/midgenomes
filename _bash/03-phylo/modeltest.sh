#!/bin/bash


#'
#' Model selection using ModelTest-NG.
#'
#'
#' Outputs:
#' - chir_modeltest
#'



. /app/.bashrc
conda activate phylo-env

# FYI, it's not much more efficient to use >18 threads:
export THREADS=$(count_threads)


export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_modeltest
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz


export ALIGNS_PART=trim_aligns.partition
cp ${TARGET}/${ALIGNS_PART}.gz ./ && gunzip ${ALIGNS_PART}.gz

# Change formatting to match expected input:
sed -i 's/WAG,/AA,/g' ${ALIGNS_PART}


#'
#' -t ml
#' Optimize topology and branch lengths for each model
#'
#' -d aa
#' amino acid
#'
#' -o modeltest-ng.log
#' pipes the output into a file
#'
#' -q ${ALIGNS_PART}
#' sets partitions file
#'
#' -p ${THREADS}
#' number of processes to use (shared memory)
#'
#' -r INTEGER
#' sets the seed for the random number generator
#'
#' -T raxml
#' sets candidate models according to a specified tool
#'

modeltest-ng -i ${CONCAT_ALIGNS} \
    -t ml \
    -d aa \
    -o modeltest-ng.log \
    -q ${ALIGNS_PART} \
    -p ${THREADS} \
    -r 832741933 \
    -T raxml

# This is stored elsewhere:
rm ${CONCAT_ALIGNS}


cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

