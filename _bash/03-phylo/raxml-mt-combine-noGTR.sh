#!/bin/bash


#'
#' Create ML tree using RAxML-NG, using ModelTest-NG models as inputs.
#'
#' Outputs:
#' - chir_raxml_mt_noGTR.tar.gz
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_mt
export PREFIX=chir_phy_mt
export ML_TREE_NAME=${PREFIX}.tree
mkdir ${OUT_DIR}
cd ${OUT_DIR}




for i in {0..9}; do
    cp ${TARGET}/chir_raxml_mt_${i}.tar.gz ./ \
        && tar -xzf chir_raxml_mt_${i}.tar.gz \
        && rm chir_raxml_mt_${i}.tar.gz \
        && cd ./chir_raxml_mt_${i} \
        && cat chir_phy_mt_${i}.raxml.mlTrees >> ../${PREFIX}.raxml.mlTrees \
        && LL=$(grep "Final LogLikelihood: " chir_phy_mt_${i}.raxml.log \
            | sed 's/Final LogLikelihood: //g') \
        && echo -e "${i}\t${LL}" >> ../${PREFIX}.raxml.logLik \
        && cd ..
    check_exit_status "copy raxml output $i" $?
done


#' How consistent are final trees?
mkdir rfdist
cd rfdist
raxml-ng --rfdist --tree ../${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF
cd ..
#' Very consistent:
#'
#' Average absolute RF distance in this tree set: 0.000000
#' Average relative RF distance in this tree set: 0.000000
#' Number of unique topologies in this tree set: 1
#'


# Run with the highest log likelihood:
best_run=$(R --vanilla --slave << EOF
llf = "${PREFIX}.raxml.logLik"
lld = read.table(llf)
cat(lld[,1][lld[,2] == max(lld[,2])])
EOF
)

cp chir_raxml_mt_${best_run}/chir_phy_mt_${best_run}.raxml.bestTree ${PREFIX}.raxml.bestTree
cp ${PREFIX}.raxml.bestTree ${ML_TREE_NAME}
mv ${ML_TREE_NAME} ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}


cd ${TARGET}
rm chir_raxml_mt_?.tar.gz


exit 0
