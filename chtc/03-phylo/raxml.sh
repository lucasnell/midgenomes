#!/bin/bash


#'
#' Phylogeny creation using RAxML-NG.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env


export OUT_DIR=chir_raxml
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export ALIGNS=cat_aligns.fasta
cp /staging/lnell/phylo/${ALIGNS}.gz ./ && gunzip ${ALIGNS}


raxml-ng --msa ${ALIGNS} --prefix chir_phy --threads ${THREADS} \
    --outgroup Anopheles_stephensi \
    --data-type AA \
    --seed 453418559 \
    --model PROTGTR+I+G \
    --bs-trees 100 \
    1> >(tee -a chir_phy.log)


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/

rm -r ${OUT_DIR}
