#!/bin/bash


#'
#' Create ML tree using IQ-TREE.
#'
#'
#' Outputs:
#' - chir_iqtree
#' - chir_ml.tree
#'



. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_iqtree
export ML_TREE_NAME=chir_ml.tree
export PREFIX=chir_phy
### mkdir ${OUT_DIR}
### cd ${OUT_DIR}
###
###
export CONCAT_ALIGNS=mafft_aligns_concat.faa
### cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
### export ALIGNS_PART=trim_aligns.partition
### cp ${TARGET}/${ALIGNS_PART}.gz ./ && gunzip ${ALIGNS_PART}.gz

cp ${TARGET}/${OUT_DIR}.tar ./ && tar -xf ${OUT_DIR}.tar
cd ${OUT_DIR}



# Convert partition file to nexus format:
export ALIGNS_PART_NEXUS=trim_aligns.nex
### echo -e "#nexus\nbegin sets;" > ${ALIGNS_PART_NEXUS}
### cat ${ALIGNS_PART} \
###     | sed 's/WAG,/  charset/g' \
###     | sed 's/$/;/' \
###     >> ${ALIGNS_PART_NEXUS}
### echo "end;" >> ${ALIGNS_PART_NEXUS}
### rm ${ALIGNS_PART}



#'
#' -m MFP+MERGE
#' Extended model selection followed by tree inference.
#' Also merge partitions to increase model fit.
#'
#' -T AUTO
#' Automatically choose number of threads to use.
#'
#' --threads-max ${THREADS}
#' Choose at max `${THREADS}` threads.
#'
#' --seqtype AA
#' Sequence type is amino acids.
#'
#' -o Mdomes
#' Outgroup is Mdomes
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
#' -B 1000
#' Run 1000 ultrafast bootstrap replicates.
#'
#' -bnni
#' With this option UFBoot will further optimize each bootstrap tree using a
#' hill-climbing nearest neighbor interchange (NNI) search based directly on
#' the corresponding bootstrap alignment.
#'
#' --alrt 1000
#' SH-like approximate likelihood ratio test
#' ([Guindon et al., 2010](https://doi.org/10.1093/sysbio/syq010))
#' to test support for individual branches
#'
#' -p FILE
#' Specify partition file for edge-proportional partition model,
#' where partitions have different evolutionary speeds.
#'






iqtree -s ${CONCAT_ALIGNS} \
    -m MFP+MERGE \
    -T AUTO \
    --threads-max ${THREADS} \
    --seqtype AA \
    -o Mdomes \
    --prefix ${PREFIX} \
    --msub nuclear \
    --seed 972232353 \
    --rclusterf 10 \
    -B 1000 \
    --bnni \
    --alrt 1000 \
    -p ${ALIGNS_PART_NEXUS} \
    | tee iqtree.log


### # This is stored elsewhere:
### rm ${CONCAT_ALIGNS}
###
### # Move ML tree to final name:
### mv ${PREFIX}.treefile ../${ML_TREE_NAME}

cd ..




tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

# rename old one in case this one fails:
mv ${TARGET}/${OUT_DIR}.tar ${TARGET}/${OUT_DIR}_old.tar

mv ${OUT_DIR}.tar.gz ${TARGET}/
### mv ${ML_TREE_NAME} ${TARGET}/

rm -r ${OUT_DIR}

