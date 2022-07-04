#!/bin/bash

#'
#' Use BUSCO to find genes in OrthoDB (diptera_odb10)
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

export TARGET=/staging/lnell/phylo/species-files


#' All species in the phylogeny, in alphabetical order:
export ALL_SPECIES=(Anopheles_stephensi \
    Belgica_antarctica \
    Chironomus_riparius \
    Chironomus_tentans \
    Chironomus_tepperi \
    Clunio_marinus \
    Culicoides_sonorensis \
    Parochlus_steinenii \
    Polypedilum_pembai \
    Polypedilum_vanderplanki \
    Propsilocerus_akamusi \
    Simulium_vittatum \
    Tanytarsus_gracilentus)


conda activate busco-env

busco --download diptera_odb10


for SPECIES in ${ALL_SPECIES[@]}; do

    export GENOME=${SPECIES}.fasta
    export ODB_DIR=${SPECIES}_odb
    export BUSCO_DIR=${SPECIES}_busco

    cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
    check_exit_status "cp genome" $?

    busco -m genome -i ${GENOME} -o ${BUSCO_DIR} \
        --cpu ${THREADS} -l $(pwd)/busco_downloads/lineages/diptera_odb10

    cp -r ${BUSCO_DIR}/run_diptera_odb10/busco_sequences/single_copy_busco_sequences ./
    mv single_copy_busco_sequences ${ODB_DIR}

    tar -czf ${ODB_DIR}.tar.gz ${ODB_DIR}
    mv ${ODB_DIR}.tar.gz ${TARGET}/
    tar -czf ${BUSCO_DIR}.tar.gz ${BUSCO_DIR}
    mv ${BUSCO_DIR}.tar.gz ${TARGET}/

    rm -r ${GENOME} ${ODB_DIR} ${BUSCO_DIR}

done


rm -r busco_downloads
conda deactivate



