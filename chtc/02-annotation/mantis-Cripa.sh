#!/bin/bash


#'
#' Functional annotations using mantis.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Memory available:
export MEMORY=$(grep "^Memory = " $_CONDOR_MACHINE_AD | sed 's/Memory\ =\ //')
# In GB, with 5 GB overhead:
MEMORY=$(( MEMORY / 1000 - 5 ))

export OUT_DIR=Cripa_mantis

mkdir working
cd working

# For databases:
mkdir dbs
# For references:
mkdir refs


eval "$(conda shell.bash hook)"
conda activate annotate-env

cp /staging/lnell/annotation/Cripa_proteins.fasta.gz ./ \
    && gunzip Cripa_proteins.fasta.gz


#' Setup MANTIS.cfg based on default from
#' https://github.com/PedroMTQ/mantis/blob/a59937d23979372c65188cf2358b69485099dc68/config/MANTIS.cfg
#' downloaded on 14 June 2022.
#' The only differences are that comments are removed and
#' default_ref_folder and resources_folder are set.
cat << EOF > MANTIS.cfg
default_ref_folder=$(pwd)/dbs
resources_folder=$(pwd)/refs
nog_ref=dmnd
nog_weight=0.8
pfam_weight=0.9
EOF


#' Download and setup databases:
mantis setup -mc MANTIS.cfg


#' Run mantis:
mantis run -mc MANTIS.cfg -i Cripa_proteins.fasta -o ${OUT_DIR} -od 315571 \
    --verbose_kegg_matrix \
    --cores ${THREADS} --memory ${MEMORY} --hmmer_threads 2

rm -rf refs dbs


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/annotation/

cd ..
rm -r working






