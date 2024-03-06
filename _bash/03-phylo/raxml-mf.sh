#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelFinder scheme/models as inputs.
#'
#' Outputs:
#' - chir_raxml_mf
#' - chir_ml_mf.tree
#'

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mf
export ML_TREE_NAME=chir_ml_mf.tree
export PREFIX=chir_phy_mf
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=chir_modfind_raxml.part
cp ${TARGET}/${ALIGNS_PART} ./

#' RAxML requires "+F" instead of just "F" (e.g., "LG+F" instead of "LGF"),
#' but ModelFinder doesn't include the "+".
#' From manually checking the `$ALIGNS_PART` file, I know that this
#' use find-replace works for this file. This is not a universal fix.
sed -i 's/F,/+F,/g' ${ALIGNS_PART}


#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --threads ${THREADS} --extra thread-pin \
    --outgroup Mdomes \
    --data-type AA \
    --model ${ALIGNS_PART} \
    --brlen scaled \
    --seed 925483383 \
    1> >(tee -a ${PREFIX}.stdout)


# #' How consistent are final trees?
# mkdir rfdist
# cd rfdist
# raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF
# cd ..

# These are stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}

#' This is the "Best-scoring ML tree (strictly bifurcating)" (RAxML-NG docs)
#' that we'll want access to later.
if [[ -f ${PREFIX}.raxml.bestTree ]]; then
    cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME}
    mv ${ML_TREE_NAME} ${TARGET}/
else
    echo ${PREFIX}.raxml.bestTree "not found!" 1>&2
fi

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



