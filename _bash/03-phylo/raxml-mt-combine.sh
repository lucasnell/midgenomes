#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelTest-NG models as inputs.
#'
#' Outputs:
#' - chir_raxml_mt.tar.gz
#' - chir_ml.tree
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mt
export PREFIX=chir_phy_mt
export ML_TREE_NAME=chir_ml.tree
mkdir ${OUT_DIR}
cd ${OUT_DIR}


# tsv of AIC by run:
echo -e "run\taic" > ${PREFIX}.raxml.aic

for i in {0..9}; do
    cp ${TARGET}/chir_raxml_mt_${i}.tar.gz ./ \
        && tar -xzf chir_raxml_mt_${i}.tar.gz \
        && rm chir_raxml_mt_${i}.tar.gz \
        && cd ./chir_raxml_mt_${i} \
        && cat chir_phy_mt_${i}.raxml.mlTrees >> ../${PREFIX}.raxml.mlTrees \
        && aic=$(grep "AIC" chir_phy_mt_${i}.raxml.log \
            | sed 's/AIC score: //g; s/[[:space:]].*$//') \
        && echo -e "${i}\t${aic}" >> ../${PREFIX}.raxml.aic \
        && cd ..
    check_exit_status "copy raxml output $i" $?
done


#' How consistent are final trees?
mkdir rfdist
cd rfdist
raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF \
    | tee ../${PREFIX}.raxml.rfdist
check_exit_status "raxml rfdist" $?
cd ..
#' Very consistent:
#'
#' Average absolute RF distance in this tree set: 0.000000
#' Average relative RF distance in this tree set: 0.000000
#' Number of unique topologies in this tree set: 1
#'


# Run with the lowest AIC:
best_run=$(R --vanilla --slave << EOF
aic_f = "${PREFIX}.raxml.aic"
aic_data = read.table(aic_f, header = TRUE)
run = aic_data[,"run"]
aic = aic_data[,"aic"]
cat(run[aic == min(aic)])
EOF
)
check_exit_status "get best run" $?

cp chir_raxml_mt_${best_run}/chir_phy_mt_${best_run}.raxml.bestTree ${PREFIX}.raxml.bestTree \
    && cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME} \
    && mv ${ML_TREE_NAME} ${TARGET}/
check_exit_status "renaming, moving trees" $?

cd .. \
    && tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR} \
    && mv ${OUT_DIR}.tar.gz ${TARGET}/ \
    && rm -r ${OUT_DIR}
check_exit_status "compressing, moving output dir" $?

# cd ${TARGET}
# rm chir_raxml_mt_?.tar.gz


exit 0
