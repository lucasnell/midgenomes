#!/bin/bash


#'
#' Sequence alignment using MAFFT.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

conda activate phylo-env


export OUT_DIR=raxml
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export ALIGNS=cat_aligns.fasta
cp /staging/lnell/phylo/${ALIGNS}.gz ./ && gunzip ${ALIGNS}


raxml-ng --msa ${ALIGNS} --prefix chir_phy --threads ${THREADS} \
    --outgroup Anopheles_stephensi \
    --data-type AA \
    --seed 453014559 \
    --model PROTGTR+I+G \
    --tree pars{25},rand{25} \
    1> >(tee -a chir_phy.out) \
    2> >(tee -a chir_phy.err >&2)





# LEFT OFF:
# WARNING: Sequences Polypedilum_vanderplanki and Propsilocerus_akamusi are exactly identical!

