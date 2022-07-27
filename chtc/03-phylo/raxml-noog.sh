#!/bin/bash


#'
#' Create ML tree using RAxML-NG - with NO OutGroup
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_noog
export ML_TREE_NAME=chir_ml_noog.tree
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz



raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} --threads ${THREADS} \
    --data-type AA \
    --model LG+I+G \
    --seed 1049211147 \
    1> >(tee -a ${PREFIX}.stdout)
#' With 32 threads on CHTC: Elapsed time: 10646.328 seconds (~2 hr 57 min)


#' How consistent are final trees?
mkdir rfdist
cd rfdist
raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF
#' Doing this interactively found only 1 unique topology, so I'll continue.
cd ..


rm ${CONCAT_ALIGNS}

#' This is the "Best-scoring ML tree (strictly bifurcating)" (RAxML-NG docs)
#' that we'll want access to later.
cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME}
mv ${ML_TREE_NAME} ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



