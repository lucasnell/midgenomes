#!/bin/bash


#'
#' Use IQ-TREE's ModelFinder to merge partitions.
#'
#'
#' Outputs:
#' - chir_merge
#' - chir_merge.partition
#'

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo


export OUT_DIR=chir_merge
export MF_PREFIX=chir_merge
export PART_FILE=${OUT_DIR}.partition

mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
export ALIGNS_PART=trim_aligns.partition
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ \
    && gunzip ${CONCAT_ALIGNS}.gz \
    && cp ${TARGET}/${ALIGNS_PART}.gz ./ \
    && gunzip ${ALIGNS_PART}.gz
check_exit_status "copy, gunzip aligns and partitions" $?

# Change 'AA' to 'TESTMERGEONLY'
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
#' --mset LG
#' restrict to just the LG model bc we only want to merge partitions.
#'
#' --mrate G
#' List of rate heterogeneity among sites. This specifies discrete Gamma rates.
#'
#' --seed INTEGER
#' Set seed for reproducible runs.
#'
#' --rclusterf 10
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
    --mset LG \
    --mrate G \
    --seed 1531525538 \
    --rclusterf 10 \
    -p ${ALIGNS_PART} \
    > ${MF_PREFIX}.stdout


# These stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}


# Move partition file to final name:
cp ${MF_PREFIX}.best_scheme ../${PART_FILE}
# Structure of this output file:
# LGF, 10441at7147_11408at7147 = 12841-13612, 36487-38168

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${PART_FILE} ${TARGET}/

rm -r ${OUT_DIR}


