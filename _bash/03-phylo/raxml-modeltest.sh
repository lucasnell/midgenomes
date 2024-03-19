#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelTest-NG models as inputs.
#'
#' Outputs:
#' - chir_raxml_mt
#' - chir_ml_mt.tree
#'

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mt
export ML_TREE_NAME=chir_ml_mt.tree
export PREFIX=chir_phy_mt
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=chir_modeltest.partition
cp ${TARGET}/${ALIGNS_PART} ./

#' This is important bc RAxML-NG assumes model 'GTR' refers to the
#' DNA version, but we're using the protein version:
sed -i 's/^GTR/PROTGTR/g' ${ALIGNS_PART}

#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --threads ${THREADS} \
    --extra thread-pin \
    --outgroup Mdomes \
    --data-type aa \
    --model ${ALIGNS_PART} \
    --brlen scaled \
    --seed 414112503 \
    1> >(tee -a ${PREFIX}.stdout)
status=$?

#' Sometimes we still get the 'CPU core oversubscription' error by chance,
#' so if that happens, we'll try it again with `--extra thread-nopin` to see
#' if that fixes the problem.
if [ "$status" != "0" ]; then
    if grep -Fq 'ERROR: CPU core oversubscription detected!' ${PREFIX}.raxml.log; then
        mv ${CONCAT_ALIGNS} ${ALIGNS_PART} ../ \
            && rm * \
            && mv ../${CONCAT_ALIGNS} ../${ALIGNS_PART} ./
        raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
            --threads ${THREADS} \
            --extra thread-nopin \
            --outgroup Mdomes \
            --data-type AA \
            --model ${ALIGNS_PART} \
            --seed 414112503 \
            1> >(tee -a ${PREFIX}.stdout)
        check_exit_status "raxml-ng thread-nopin" $?
    else
        check_exit_status "raxml-ng" $status
    fi
fi



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



