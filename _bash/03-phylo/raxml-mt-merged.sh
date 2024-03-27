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

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mt_merged
export ML_TREE_NAME=chir_ml_merged.tree
export PREFIX=chir_phy_mt_merged
mkdir ${OUT_DIR}
cd ${OUT_DIR}

# File in TARGET that will allow me to observe progress:
export PROG_FILE=${TARGET}/raxml-mt-merged.txt


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=chir_modeltest_merged.partition
cp ${TARGET}/${ALIGNS_PART} ./


#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
run_raxml () {
    raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
        --threads ${THREADS} \
        --extra ${1} \
        --outgroup Mdomes \
        --data-type aa \
        --model ${ALIGNS_PART} \
        --brlen scaled \
        --seed 11843774 \
        | tee -a ${PREFIX}.stdout
    local stts=$?
    return $stts
}
export THREAD_PIN="thread-pin"
run_raxml ${THREAD_PIN}
export status=$?


#' Sometimes we still get the 'CPU core oversubscription' error by chance,
#' so if that happens, we'll try it again with `--extra thread-nopin` to see
#' if that fixes the problem.
if [ "$status" != "0" ]; then
    if grep -Fq 'ERROR: CPU core oversubscription detected!' ${PREFIX}.raxml.log; then
        mv ${CONCAT_ALIGNS} ${ALIGNS_PART} ../ \
            && rm * \
            && mv ../${CONCAT_ALIGNS} ../${ALIGNS_PART} ./
        export THREAD_PIN="thread-nopin"
        run_raxml ${THREAD_PIN}
        check_exit_status "raxml-ng thread-nopin" $?
    else
        check_exit_status "raxml-ng" $status
    fi
else
    check_exit_status "raxml-ng" $status
fi


cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME}
mv ${ML_TREE_NAME} ${TARGET}/


#' How consistent are final trees?
mkdir rfdist
cd rfdist
raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF \
    | tee -a ${PROG_FILE} ../${PREFIX}.raxml.rfdist
check_exit_status "raxml rfdist" $?
cd ..



echo -e "\n\nStarted bootstrapping @ " $(date) >> ${PROG_FILE}


export BOOT_PREFIX=chir_phy_boot
export BOOT_DIR=chir_phy_boot

mkdir ${BOOT_DIR}
cd ${BOOT_DIR}

#' I'm making it perform the MRE-based bootstopping test every 50 replicates
#' and terminate once it surpasses the cutoff value of 0.01 (a strict value),
#' for a maximum of 100 trees.
#'
raxml-ng --bootstrap --msa ../${CONCAT_ALIGNS} --prefix ${BOOT_PREFIX} \
    --bs-trees autoMRE{100} \
    --threads ${THREADS} \
    --extra ${THREAD_PIN} \
    --bs-cutoff 0.01 \
    --outgroup Mdomes \
    --data-type aa \
    --model ../${ALIGNS_PART} \
    --brlen scaled \
    --seed 158939107 \
    | tee ${BOOT_PREFIX}.stdout
check_exit_status "raxml bootstrap" $?


echo -e "\n\nFinished bootstrapping @ " $(date) >> ${PROG_FILE}



#' Now look at branch support
raxml-ng --support \
    --tree ../${PREFIX}.raxml.bestTree \
    --bs-trees ${BOOT_PREFIX}.raxml.bootstraps \
    --outgroup Mdomes \
    --prefix chir_supp --threads ${THREADS} --bs-metric fbp,tbe \
    | tee chir_raxml_support.stdout

cd ..


# These are stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


