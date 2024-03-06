#!/bin/bash


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DEPRECATED - NOT NECESSARY TO RUN THIS GIVEN THAT PARTITIONING BY GENE WORKS
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#'
#' Create ML tree using RAxML-NG, using no partitioning.
#'
#' Requires an integer argument indicating what type of model to use:
#' 0 - LG
#' 1 - Q.insect
#' 2 - WAG
#'
#' Outputs:
#' - chir_raxml_simp_${RAXML_MODEL}
#' - chir_ml_simp_${RAXML_MODEL}.tree
#'
#' Usage:
#' raxml.sh RAXML_MODEL_INDEX
#'

if (( $# != 1 )); then
    echo "ERROR: Exactly one argument required for raxml.sh" 1>&2
    exit 1
fi

RAXML_MODEL_INDEX=$1

if ! [[ $RAXML_MODEL_INDEX =~ ^[0-9]+$ ]] || (( RAXML_MODEL_INDEX < 0 )) ||
        (( RAXML_MODEL_INDEX >= 3 )); then
    echo -n "ERROR: First argument should be an integer 0-2. " 1>&2
    echo "Yours is '${RAXML_MODEL_INDEX}'." 1>&2
    exit 1
fi



ALL_RAXML_MODELS=("LG" "Q.insect" "WAG")
export RAXML_MODEL=${ALL_RAXML_MODELS[$RAXML_MODEL_INDEX]}

ALL_RAXML_SEEDS=(1516919124 304060730 1385024301)
export RAXML_SEED=${ALL_RAXML_SEEDS[$RAXML_MODEL_INDEX]}

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_simp_${RAXML_MODEL}
export ML_TREE_NAME=chir_ml_simp_${RAXML_MODEL}.tree
export PREFIX=chir_phy_simp_${RAXML_MODEL}
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz


#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --threads ${THREADS} --extra thread-pin \
    --outgroup Mdomes \
    --data-type AA \
    --model $(echo ${RAXML_MODEL}"+I+G") \
    --seed ${RAXML_SEED} \
    1> >(tee -a ${PREFIX}.stdout)
check_exit_status "raxml-ng" $?
#' With 32 threads on CHTC: Elapsed time: 31172.097 seconds (~8 hr 40 min)


# #' How consistent are final trees?
# mkdir rfdist
# cd rfdist
# raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF
# cd ..


rm ${CONCAT_ALIGNS}

#' This is the "Best-scoring ML tree (strictly bifurcating)" (RAxML-NG docs)
#' that we'll want access to later.
cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME}
mv ${ML_TREE_NAME} ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



