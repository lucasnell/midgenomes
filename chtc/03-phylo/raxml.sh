#!/bin/bash


#'
#' Create ML tree using RAxML-NG.
#'




export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml
export ML_TREE_NAME=chir_ml.tree
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz


# Started 10:55am
date > start-end-times.out
#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --threads ${THREADS} --extra thread-pin \
    --outgroup Mdomes \
    --data-type AA \
    --model LG+I+G \
    --seed 1171225701 \
    1> >(tee -a ${PREFIX}.stdout)
date >> start-end-times.out
#' With 16 threads on CHTC: Elapsed time: 18354.721 seconds (~5 hr 6 min)
#' With 32 threads on CHTC: Elapsed time:  seconds (~)


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



