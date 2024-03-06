#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using partitioning from ModelFinder, but
#' manually fitting substition model.
#'
#' Requires an integer argument indicating what type of model to use:
#' 0 - LG
#' 1 - JTT
#' 2 - WAG
#' 3 - Q.insect
#' 4 - Q.pfam
#'
#' Outputs:
#' - chir_raxml_mfmerge_${RAXML_MODEL}
#' - chir_ml_mfmerge_${RAXML_MODEL}.tree
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
        (( RAXML_MODEL_INDEX > 4 )); then
    echo -n "ERROR: First argument should be an integer 0-4. " 1>&2
    echo "Yours is '${RAXML_MODEL_INDEX}'." 1>&2
    # exit 1
fi


ALL_RAXML_MODELS=("LG" "JTT" "WAG" "Q.insect" "Q.pfam")
export RAXML_MODEL=${ALL_RAXML_MODELS[$RAXML_MODEL_INDEX]}

ALL_RAXML_SEEDS=(407078463 314333995 851292978 1725688976 958824602)
export RAXML_SEED=${ALL_RAXML_SEEDS[$RAXML_MODEL_INDEX]}

. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mfmerge_${RAXML_MODEL}
export ML_TREE_NAME=chir_ml_mfmerge_${RAXML_MODEL}.tree
export PREFIX=chir_phy_mfmerge_${RAXML_MODEL}
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=chir_modfind_raxml.part
cp ${TARGET}/${ALIGNS_PART} ./

# Change models to the one we're choosing:
FULL_MODEL=$(echo ${RAXML_MODEL}"+I+G")
sed -i "s/^[^,]*,/${FULL_MODEL},/" ${ALIGNS_PART}



#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --search --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --threads ${THREADS} \
    --extra thread-pin \
    --outgroup Mdomes \
    --data-type AA \
    --model ${ALIGNS_PART} \
    --seed ${RAXML_SEED} \
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
            --seed ${RAXML_SEED} \
            1> >(tee -a ${PREFIX}.stdout)
        check_exit_status "raxml-ng thread-nopin" $?
    else
        check_exit_status "raxml-ng" $status
    fi
fi


# Extract AICc:
aicc=$(grep "AIC score:" ${PREFIX}.raxml.log | sed 's/^[^\/]*\/ //; s/ \/.*//; s/AICc score: //g')
# Write to file for simpler comparison:
echo -e "${FULL_MODEL}\t${aicc}" >> ${TARGET}/${OUT_DIR/_${RAXML_MODEL}/}_AICc.tsv


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



