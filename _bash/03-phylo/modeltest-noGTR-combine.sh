#!/bin/bash


#'
#' Combine model selection using ModelTest-NG.
#'
#'
#' Outputs:
#' - chir_modeltest_noGTR.tar.gz
#' - chir_modeltest_noGTR.partition
#'

. /app/.bashrc

export TARGET=/staging/lnell/phylo

export PREFIX=chir_modeltest_noGTR
export OUT_DIR=${PREFIX}
export OUT_PART=${PREFIX}.partition
mkdir ${OUT_DIR}
cd ${OUT_DIR}


for i in {0..9}; do
    cp ${TARGET}/${PREFIX}_${i}.tar.gz ./
    check_exit_status "null" $?
    tar -xzf ${PREFIX}_${i}.tar.gz
    rm ${PREFIX}_${i}.tar.gz
done

#' rename files to remove unnecessary '.log' and combine the files
#' modeltest_${i}.part.aic into one
for i in {0..9}; do
    cd ${PREFIX}_${i}
    sed -i "s/modeltest_${i}.log/modeltest_${i}/g" modeltest_${i}.log.log
    sed -i "s/modeltest_${i}.log/modeltest_${i}/g" modeltest_${i}.log.out
    for f in modeltest_${i}.log*; do
        mv $f ${f/.log/}
    done
    cat modeltest_${i}.part.aic >> ../${OUT_PART}
    cd ..
done

cp ${OUT_PART} ${TARGET}/

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

