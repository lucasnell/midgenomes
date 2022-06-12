#!/bin/bash

#'
#' Use BUSCO to find genes in OrthoDB (diptera_odb10)
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


export TARGET=/staging/lnell/phylo



cd /staging/lnell/phylo
export ALL_SPECIES=($(ls *.fasta.gz | sed 's/.fasta.gz//g'))
cd $_CONDOR_SCRATCH_DIR


conda activate busco-env

busco --download diptera_odb10


for SPECIES in ${ALL_SPECIES[@]}; do

    export GENOME=${SPECIES}.fasta
    export OUT_DIR=${SPECIES}_odb

    cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
    check_exit_status "cp genome" $?

    busco -m genome -i ${GENOME} -o busco \
        --cpu ${THREADS} -l $(pwd)/busco_downloads/lineages/diptera_odb10

    mv busco/run_diptera_odb10/busco_sequences/single_copy_busco_sequences ./
    mv single_copy_busco_sequences ${OUT_DIR}

    tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ${OUT_DIR}.tar.gz ${TARGET}/
    rm -r ${GENOME} ${OUT_DIR} busco

done


rm -r busco_downloads
conda deactivate



