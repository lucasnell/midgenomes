#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelTest-NG models as inputs.
#' This also does bootstrapping and computes branch support.
#'
#'
#' Outputs:
#' - chir_raxml_mt_merged.tar.gz
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)
export WORKERS=$(( THREADS / 4 ))

export TARGET=/staging/lnell/phylo


export PREFIX=chir_raxml_boot
export OUT_DIR=chir_raxml_mt_merged_boot

mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
export ALIGNS_PART=chir_modeltest_merged.partition
export ML_TREE_NAME=chir_ml_merged.tree

cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
cp ${TARGET}/${ALIGNS_PART} ./
cp ${TARGET}/${ML_TREE_NAME} ./


#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
run_raxml_boot () {
    #'
    #' I'm making it perform the MRE-based bootstopping test every 50 replicates
    #' and terminate once it surpasses the cutoff value of 0.01 (a strict value),
    #' up to a maximum of 1000 trees.
    #'
    raxml-ng --bootstrap --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
        --bs-trees autoMRE{1000} \
        --threads ${THREADS} \
        --workers auto{} \
        --extra ${1} \
        --bs-cutoff 0.01 \
        --outgroup Mdomes \
        --data-type aa \
        --model ${ALIGNS_PART} \
        --brlen scaled \
        --seed 158939107 \
        | tee ${PREFIX}.stdout
    local stts=$?
    return $stts
}


run_raxml_boot "thread-pin"
export status=$?


#' Sometimes we still get the 'CPU core oversubscription' error by chance,
#' so if that happens, we'll try it again with `--extra thread-nopin` to see
#' if that fixes the problem.
if [ "$status" != "0" ]; then
    if grep -Fq 'ERROR: CPU core oversubscription detected!' ${PREFIX}.raxml.log; then
        mv ${CONCAT_ALIGNS} ${ALIGNS_PART} ../ \
            && rm * \
            && mv ../${CONCAT_ALIGNS} ../${ALIGNS_PART} ./
        run_raxml_boot "thread-nopin"
        check_exit_status "raxml-ng boot thread-nopin" $?
    else
        check_exit_status "raxml-ng boot" $status
    fi
else
    check_exit_status "raxml-ng boot" $status
fi





#' Now look at branch support
mkdir chir_supp
cd chir_supp
raxml-ng --support \
    --tree ../${ML_TREE_NAME} \
    --bs-trees ../${PREFIX}.raxml.bootstraps \
    --outgroup Mdomes \
    --prefix chir_raxml_supp \
    --threads auto{${THREADS}} \
    --bs-metric fbp,tbe \
    | tee chir_raxml_supp.stdout
check_exit_status "raxml-ng support" $?
cd ..


# These are stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART} ${ML_TREE_NAME}


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


