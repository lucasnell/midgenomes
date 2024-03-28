#!/bin/bash


#'
#' Find orthogroups using OrthoFinder.
#' Ran this using interactive job with 48 threads, 64 GB RAM, 150 GB disk.
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_orthofinder
mkdir ${OUT_DIR}
cd ${OUT_DIR}

#' Time tree in NEWICK format from MCMCTree output.
export SPECIES_TREE=time-tree.nwk
cp ${TARGET}/phylo/${SPECIES_TREE} ./
check_exit_status "copy time tree" $?


export PROT_FOLDER=chir_proteins

cp -r ${TARGET}/proteins ./ \
    && mv proteins ${PROT_FOLDER} \
    && cd ${PROT_FOLDER}
check_exit_status "moving, renaming protein folder" $?



# Takes just a couple of minutes
for f in *.faa.gz; do
    g=$(echo ${f%.gz} | sed 's/_proteins//g')
    spp=${f/_proteins.faa.gz/}
    echo $spp >> longest-isoform.log
    longest-isoform.py $f 1> ${g} 2>> longest-isoform.log
    check_exit_status "only longest isoforms - $spp" $?
    echo -e "\n\n" >> longest-isoform.log
    rm $f
done

mv longest-isoform.* ../
cd ..



orthofinder -f ${PROT_FOLDER} -t ${THREADS} -a $(( THREADS / 4 )) \
    -s ${SPECIES_TREE} \
    2>&1 \
    | tee orthofinder.log
check_exit_status "run OrthoFinder" $?

# Move OrthoFinder output out of the proteins folder and rename:
cd ${PROT_FOLDER} \
    && mv OrthoFinder OrthoFinder_tmp \
    && cd OrthoFinder_tmp \
    && mv Results_* orthofinder-output \
    && mv orthofinder-output ../../ \
    && cd .. \
    && rm -r OrthoFinder_tmp \
    && cd ..
check_exit_status "move, rename OrthoFinder output" $?

# Remove large working directory folder:
rm -r ${OUT_DIR}/orthofinder-output/WorkingDirectory

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR} \
    && mv ${OUT_DIR}.tar.gz ${TARGET}/
check_exit_status "move output to target directory" $?

rm -r ${OUT_DIR}

