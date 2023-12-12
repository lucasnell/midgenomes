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

export PROT_FOLDER=chir_proteins

cp -r /staging/lnell/proteins ./ \
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



#' Create simplified time tree in NEWICK format from MCMCTree output.
#' Use this here in OrthoFinder and save for potentially using elsewhere.
export SPECIES_TREE_MCMCTREE=chir_mcmctree.tre
export SPECIES_TREE=chir_mcmctree.nwk
cat /staging/lnell/phylo/${SPECIES_TREE_MCMCTREE} \
    | sed -e 's/\[[^][]*\]//g' \
    | grep "UTREE" \
    | sed 's/[^(]*//' \
    > ${SPECIES_TREE}
check_exit_status "create simple newick time tree" $?

#' This forces the NEWICK tree to be ultrametric. Taken from
#' https://github.com/PuttickMacroevolution/MCMCtreeR/blob/2330e7a9916c3929513ee217d3854be993965f6b/R/readMCMCTree.R#L53-L70
#'
R --vanilla << EOF
library(ape)
phy <- read.tree("${SPECIES_TREE}")
outer <- phy\$edge[,2]
inner <- phy\$edge[,1]
totalPath <- c()
for(i in which(outer<=Ntip(phy))) {
    start <- i
    end <- inner[start]
    edgeTimes <- phy\$edge.length[start]
    while(end != inner[1]) {
        start <- which(outer == end)
        end <- inner[start]
        edgeTimes <- c(edgeTimes, phy\$edge.length[start])
    }
    totalPath <- c(totalPath, sum(edgeTimes))
}
addLength <- max(totalPath) - totalPath
phy\$edge.length[which(outer <= Ntip(phy))] <- phy\$edge.length[
    which(outer <= Ntip(phy))] + addLength
write.tree(phy,"${SPECIES_TREE}")
EOF
check_exit_status "make simple newick time tree ultrametric" $?

cp ${SPECIES_TREE} /staging/lnell/phylo/

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

