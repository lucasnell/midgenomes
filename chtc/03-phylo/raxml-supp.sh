#!/bin/bash


#'
#' Use bootstrapped trees to evaluate branch support using RAxML-NG.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_supp
export ML_TREE_NAME=chir_ml.tree
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp ${TARGET}/${ML_TREE_NAME}.gz ./ && gunzip ${ML_TREE_NAME}

#' combine bootstrapped trees:
export ALL_BOOTS=chir_phy_all.raxml.bootstraps
for t in ${TARGET}/chir_raxml_boot_*; do
    tar -xzf ${t} -C ./
    boot_dir=$(basename ${t} | sed 's/\.tar\.gz//g')
    cat ${boot_dir}/chir_phy.raxml.bootstraps >> ${ALL_BOOTS}
    rm -r ${boot_dir}
    unset boot_dir
done


#' Did the boostraps converge? (--bs-cutoff 0.01 is extra strict)
raxml-ng --bsconverge --bs-trees ${ALL_BOOTS} --prefix chir_bsconverge \
    --seed 387034756 --threads ${THREADS} --bs-cutoff 0.01



#' If that looks good, look at branch support
raxml-ng --support --tree ${ML_TREE_NAME} --bs-trees ${ALL_BOOTS} \
    --outgroup Anopheles_stephensi \
    --prefix chir_supp --threads ${THREADS} --bs-metric fbp,tbe


rm ${ML_TREE_NAME}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
