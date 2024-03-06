#!/bin/bash


#'
#' Create ML tree using IQ-TREE, using a non-reversible model to infer
#' a *rooted* tree and using Mdomes as outgroup.
#'
#'
#' Outputs:
#' - chir_iqtree_nonrev_o
#' - chir_ml_nonrev_o.tree
#'



. /app/.bashrc
conda activate phylo-env

# FYI, it's not much more efficient to use >18 threads:
export THREADS=$(count_threads)


export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_iqtree_nonrev_o
export ML_TREE_NAME=chir_ml_nonrev_o.tree
export PREFIX=chir_phy_nonrev_o
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz

# From a previous run using ModelFinder to merge partitions:
export ALIGNS_PART_NEXUS=chir_modfind_parts.nex
cp ${TARGET}/${ALIGNS_PART_NEXUS} ./


#' -m MFP
#' Extended model selection followed by tree inference.
#' NO merging of partitions.
#'
#' --mset NONREV,NQ.insect,NQ.pfam
#' In model selection, choose from these time non-reversible amino acid models:
#' the most general amino-acid model (NONREV),
#' an empirically estimated time model for insects (NQ.insect), or
#' an empirically estimated time model for many taxa (NQ.pfam).
#'
#' -o Mdomes
#' Outgroup is Mdomes
#'
#' -T ${THREADS}
#' Use `${THREADS}` threads.
#'
#' --prefix ${PREFIX}
#' Prefix for output files is `${PREFIX}`
#'
#' --seed INTEGER
#' Set seed for reproducible runs.
#'
#' -B 1000
#' Run 1000 ultrafast bootstrap replicates.
#'
#' --bnni
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
#' --test 10000 --test-au
#' Re-root the tree on every branch and perform tree topology tests to compare
#' the log-likelihoods of the trees being rooted on every branch of the ML tree.
#' Also perform several tree topology tests including the
#' approximately-unbiased (AU) test for the output tree.
#'



iqtree -s ${CONCAT_ALIGNS} \
    -m MFP \
    --mset NONREV,NQ.insect,NQ.pfam \
    -o Mdomes \
    -T ${THREADS} \
    --prefix ${PREFIX} \
    --seed 2036076548 \
    -B 1000 \
    --bnni \
    --alrt 1000 \
    -p ${ALIGNS_PART_NEXUS} \
    --test 10000 --test-au \
    | tee iqtree.log




# This is stored elsewhere:
rm ${CONCAT_ALIGNS}

# Move ML tree to final name:
cp ${PREFIX}.treefile ../${ML_TREE_NAME}

cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/
mv ${ML_TREE_NAME} ${TARGET}/

rm -r ${OUT_DIR}

