#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelTest-NG models as inputs.
#' One run of this script performs tree searches using just two starting trees,
#' 1 random and 1 parsimony-based.
#' You'll have to run this 10 times (with different seeds) to replicate
#' the default in RAxML-NG of 10 random and 10 parsimony-based starting trees.
#'
#' Requires an integer argument for which rep this is (0-9).
#'
#' Outputs:
#' - chir_raxml_mt_noGTR_${index}.tar.gz
#'


if (( $# != 1 )); then
    echo "ERROR: Exactly one argument required for raxml-mt-sep-noGTR.sh" 1>&2
    exit 1
fi

export index=$1

if ! [[ $index =~ ^[0-9]+$ ]] || (( index < 0 )) ||
        (( index > 9 )); then
    echo -n "ERROR: First argument should be an integer 0-9. " 1>&2
    echo "Yours is '${index}'." 1>&2
    exit 1
fi



all_seeds=(1082733730 3061676 1140808739 195385868 358982958
           87345885 2121316106 1606111671 434000947 868155269)
export rng_seed=${all_seeds[$index]}



. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mt_noGTR_${index}
export PREFIX=chir_phy_mt_${index}
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
export ALIGNS_PART=chir_modeltest_noGTR.partition
cp ${TARGET}/${ALIGNS_PART} ./


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
    --seed ${rng_seed} \
    --tree pars{1},rand{1} \
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
            --data-type aa \
            --model ${ALIGNS_PART} \
            --brlen scaled \
            --seed ${rng_seed} \
            --tree pars{1},rand{1} \
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


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



