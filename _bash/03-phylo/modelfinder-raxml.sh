#!/bin/bash


#'
#' Find optimum model and partition scheme using IQ-TREE's ModelFinder,
#' for use in RAxML-NG.
#'
#'
#' Outputs:
#' - chir_modfind_raxml
#' - chir_modfind_raxml.part
#'



. /app/.bashrc
conda activate phylo-env

# FYI, it's not much more efficient to use >18 threads:
export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_modfind_raxml
export PREFIX=${OUT_DIR}
export PART_FILE=chir_modfind_raxml.part
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=trim_aligns.partition
cp ${TARGET}/${ALIGNS_PART}.gz ./ && gunzip ${ALIGNS_PART}.gz


# Convert partition file to nexus format:
export ALIGNS_PART_NEXUS=trim_aligns.nex
echo -e "#nexus\nbegin sets;" > ${ALIGNS_PART_NEXUS}
cat ${ALIGNS_PART} \
    | sed 's/WAG,/  charset/g' \
    | sed 's/$/;/' \
    >> ${ALIGNS_PART_NEXUS}
echo "end;" >> ${ALIGNS_PART_NEXUS}
rm ${ALIGNS_PART}
unset -v ALIGNS_PART


#'
#' -m MF+MERGE
#' Extended model selection and merge partitions to increase model fit.
#' Does NOT include tree inference.
#'
#' --mset raxml
#' Restrict search to models supported by raxml.
#'
#' -T ${THREADS}
#' Use `${THREADS}` threads.
#'
#' --seqtype AA
#' Sequence type is amino acids.
#'
#' --prefix ${PREFIX}
#' Prefix for output files is `${PREFIX}`
#'
#' --msub nuclear
#' restrict to those AA models designed for nuclear source
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
    -m MF+MERGE \
    --mset raxml \
    -T ${THREADS} \
    --seqtype AA \
    --prefix ${PREFIX} \
    --msub nuclear \
    --seed 311949017 \
    --rclusterf 10 \
    -p ${ALIGNS_PART_NEXUS} \
    | tee modelfinder.log


# This is stored elsewhere:
rm ${CONCAT_ALIGNS}


# Move partition file to final name:
cp ${PREFIX}.best_scheme ../${PART_FILE}

cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/
mv ${PART_FILE} ${TARGET}/

rm -r ${OUT_DIR}

