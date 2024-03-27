#!/bin/bash


#'
#' Ultrametric timetree based on ML tree and fossil records using mcmctree
#' in PAML.
#'
#' NOTE: PAML is mostly serial, so just `$N_MCMC_RUNS` threads is fine to
#' run the final steps in parallel.
#'

. /app/.bashrc
conda activate phylo-env

#' If you ever increase this, change `SEEDS` below.
export N_MCMC_RUNS=4

#' Ideally you provide a number of threads equal to the number of eventual
#' MCMCtree runs.
export THREADS=$(count_threads)

#' No need to use more than `N_MCMC_RUNS` threads:
if (( THREADS > N_MCMC_RUNS )); then
    THREADS=$N_MCMC_RUNS
fi

#' Make sure this array is at least as long as the number of replicate
#' MCMCtree runs + 1.
export SEEDS=(264611523 791063585 720616543 2049853014 274801972 218694456
              827775792 1408775991 467020590 889086769)
if (( N_MCMC_RUNS+1 > ${#SEEDS[@]} )); then
    echo "ERROR: not enough seeds provided" 1>&2
    # Exit but not in an interactive session:
    if [[ $- != *i* ]]; then exit 1; fi
fi

export TARGET=/staging/lnell/phylo

mkdir working
cd working

export OUT_DIR=chir_mcmctree
export OUT_TREE=chir_mcmctree.nwk
export OUT_PREFIX=chir_mcmctree

mkdir ${OUT_DIR}
cd ${OUT_DIR}



export HESSIANS_DIR=chir_mcmctree_hessians
tar -xzf ${TARGET}/${HESSIANS_DIR}.tar.gz -C ./

export ALIGN_PHYLIP_PART=mafft_aligns_part.phylip
export MCMC_INPUT_TREE=mcmctree_in_tree.nwk

mv ${HESSIANS_DIR}/${ALIGN_PHYLIP_PART} ./
mv ${HESSIANS_DIR}/${MCMC_INPUT_TREE} ./
mv ${HESSIANS_DIR}/hessians/mcmctree.ctl ./

export N_GENES=$(grep "^ " ${ALIGN_PHYLIP_PART} | wc -l)


for ((i=1; i<=N_GENES; i++)); do
    ii=$(printf "%04d\n" $i)
    grep "^tree length" ./${HESSIANS_DIR}/hessians/part_${ii}/tmp${ii}.out \
        | sed 's/.*=//g' \
        | xargs \
        >> tree_lengths.txt
    phylo_string=$()
R --vanilla --slave << EOF
t = ape::read.tree(text = "${phylo_string}")
mean(ape::node.depth.edgelength(t)[1:t$Nnode])
EOF
    # LEFT OFF -- DOESN'T WORK FOR ALL GENES. MAYBE TRY READING PHYLOGENY
    # DIRECTLY?
done



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# Move Hessian outputs to in.BV file:
for ((i=1; i<=N_GENES; i++)); do
    cat ./${HESSIANS_DIR}/hessians/part_$(printf "%04d\n" $i)/rst2 >> in.BV
    check_exit_status "null" $?
done


#' --------------------
#' Create directories for the different runs of MCMCtree, each with a
#' different starting seed.
#' --------------------

INPUT_DIR=$(pwd)
for ((i=1; i<=${N_MCMC_RUNS}; i++)); do
    MCMC_DIR=mcmc_${i}
    mkdir ${MCMC_DIR}
    cat mcmctree.ctl \
        | sed "s/^usedata.*/usedata\ =\ 2 in.BV/g; s/^seed.*/seed\ =\ ${SEEDS[$i]}/g" \
        | sed "s|^seqfile.*|seqfile\ =\ ${INPUT_DIR}/${ALIGN_PHYLIP_PART}|g" \
        | sed "s|^treefile.*|treefile\ =\ ${INPUT_DIR}/${MCMC_INPUT_TREE}|g" \
        | sed "s|^mcmcfile.*|mcmcfile\ =\ ${OUT_PREFIX}_${i}.txt|g" \
        | sed "s|^outfile.*|outfile\ =\ ${OUT_PREFIX}_${i}.out|g" \
        | sed "s|^burnin|burnin = 2000|g" \
        | sed "s|^sampfreq|sampfreq = 5|g" \
        | sed "s|^nsample|nsample = 50000|g" \
        > ${MCMC_DIR}/mcmctree_${i}.ctl
    cp in.BV ./${MCMC_DIR}/
    unset -v MCMC_DIR
done
unset -v INPUT_DIR



#' Simple script to run a single iteration of MCMCtree.
#' Usage:
#'  run_mcmctree.sh [RUN NUMBER]
cat << EOF > run_mcmctree.sh
#!/bin/bash
cd mcmc_\${1}
status=\$?
if (( status != 0 )); then
    exit 1
fi
mcmctree mcmctree_\${1}.ctl \\
    1> mcmctree_\${1}.stdout \\
    2> mcmctree_\${1}.stderr
status=\$?
cd ..
exit \$status
EOF
chmod +x run_mcmctree.sh


#'
#' Python code to run MCMCtree multiple times using multiple threads.
#' This step takes a while to run (around a day).
#'
python3 << EOF
import subprocess as sp
import multiprocessing as mp
import sys
import os

def work(mcmc_run):
    """Defines the work unit for one run"""
    cmd = "./run_mcmctree.sh " + str(mcmc_run)
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    n_mcmc_runs = int(os.environ["N_MCMC_RUNS"])
    tasks = [x+1 for x in range(n_mcmc_runs)]
    n_tasks = len(tasks)
    n_threads = int(os.environ["THREADS"])
    with mp.Pool(processes=n_threads) as pool:
        return_codes = pool.map(work, tasks)
        return_codes = [str(x) for x in return_codes]
        print("return codes:\n" + "\n".join(return_codes))
    sys.exit(0)

EOF



#' --------------------
#' Process output tree
#' --------------------

#' Create simplified time tree in NEWICK format from MCMCTree output.
cat ./mcmc_1/FigTree.tre \
    | sed -e 's/\[[^][]*\]//g' \
    | grep "UTREE" \
    | sed 's/[^(]*//' \
    > ${OUT_TREE}
check_exit_status "create simple newick time tree" $?

#' This forces the NEWICK tree to be ultrametric. Taken from
#' https://github.com/PuttickMacroevolution/MCMCtreeR/blob/2330e7a9916c3929513ee217d3854be993965f6b/R/readMCMCTree.R#L53-L70
#'
R --vanilla --slave << EOF
library(ape)
file_name <- Sys.getenv("OUT_TREE")
stopifnot(grepl(".nwk\$", file_name))
phy <- read.tree(file_name)
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
write.tree(phy, file_name)
EOF
check_exit_status "make simple newick time tree ultrametric" $?



#' --------------------
#' Move output files
#' --------------------

cp ${OUT_TREE} ${TARGET}/
cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working
