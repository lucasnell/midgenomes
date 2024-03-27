#!/bin/bash


#'
#' Use IQ-TREE's ModelFinder to merge partitions for use in MCMCTree.
#'
#'
#' Outputs:
#' - chir_modelfinder
#' - chir_modelfinder.part
#'

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

mkdir working
cd working

export OUT_DIR=chir_modelfinder
export MF_PREFIX=chir_mf
export PART_FILE=${OUT_DIR}.partition

mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=trim_aligns.partition
cp ${TARGET}/${ALIGNS_PART}.gz ./ && gunzip ${ALIGNS_PART}.gz

# Change 'AA' to 'MF+MERGE'
sed -i 's/AA/TESTMERGEONLY/g' ${ALIGNS_PART}

#'
#' -m TESTMERGEONLY
#' Extended model selection and merge partitions to increase model fit.
#' Does NOT include tree inference.
#'
#' -T ${THREADS}
#' Use `${THREADS}` threads.
#'
#' --seqtype AA
#' Sequence type is amino acids.
#'
#' --prefix ${MF_PREFIX}
#' Prefix for output files is `${MF_PREFIX}`
#'
#' --mset Dayhoff,DCMut,JTT,JTTDCMut,LG,mtART,mtZOA,WAG
#' restrict to the AA models present in MCMCTree
#'
#' --mrate G5
#' List of rate heterogeneity among sites.
#' This specifies a 5-category discrete Gamma
#'
#' --seed INTEGER
#' Set seed for reproducible runs.
#'
#' --rclusterf 100
#' Specify the percentage for the fast relaxed clustering
#' algorithm (Lanfear et al., 2017) to speed up the computation instead of
#' the default slow greedy algorithm. This is similar to --rcluster-percent
#' option of PartitionFinder. For example, with -rclusterf 10 only the top
#' 10% partition schemes are considered to save computations.
#'
#' -p FILE
#' Specify partition file for edge-proportional partition model,
#' where partitions have different evolutionary speeds.
#'



iqtree -s ${CONCAT_ALIGNS} \
    -m TESTMERGEONLY \
    -T ${THREADS} \
    --seqtype AA \
    --prefix ${MF_PREFIX} \
    --mset Dayhoff,DCMut,JTT,JTTDCMut,LG,mtART,mtZOA,WAG \
    --mrate G5 \
    --seed 1531525538 \
    --rclusterf 100 \
    -p ${ALIGNS_PART} \
    | tee modelfinder.log


# These stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}


# Move partition file to final name:
cp ${MF_PREFIX}.best_scheme ../${PART_FILE}

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${PART_FILE} ${TARGET}/


cd ..
rm -r working


