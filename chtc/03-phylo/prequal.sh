#!/bin/bash


#'
#' Pre-alignment quality filter with prequal
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

conda activate phylo-env

cp /staging/lnell/phylo/common_genes.tar.gz ./
tar -xzf common_genes.tar.gz
rm common_genes.tar.gz


export OUT_DIR=prequal_filt
mkdir ${OUT_DIR}


#' Script to run prequal to do filtering on a single gene.
#' For input `X.faa`, prequal automatically outputs 3 files:
#'     X.faa.filtered, X.faa.filtered.PP, and X.faa.warning.
#'     This script also adds X.faa.stdout.
#' Usage:
#'      ./one_filter.sh [*.faa FILE]
cat << EOF > one_filter.sh
#!/bin/bash
cp \$1 ./${OUT_DIR}
IN=\$(basename \$1)
cd ${OUT_DIR}
LOG=\${IN}.stdout
prequal \${IN} > \${LOG}
rm \${IN}
cd ..
EOF
chmod +x one_filter.sh


python3 << EOF
#!/usr/bin/env python3
import glob
import subprocess as sp
import multiprocessing as mp
import sys

def work(fa_file):
    """Defines the work unit on an input file"""
    cmd = "./one_filter.sh " + fa_file
    ret = sp.run(cmd, shell = True)
    return ret

tasks = glob.glob("common_genes/*.faa")

n_tasks = len(tasks)
with mp.Pool(processes=${THREADS}) as pool:
    for i, _ in enumerate(pool.imap_unordered(work, tasks), 1):
        sys.stderr.write('\rdone {0:%}'.format(i/n_tasks))
sys.stderr.write('\n')

sys.exit(0)

EOF




#'
#' Exclude any genes where one or more species had sequences fully removed
#'
cd ${OUT_DIR}
mkdir exclude
FULLY_RM=($(egrep -lir --include=*.warning "Fully removed" . \
    | sed 's/\.faa\.warning/\*/g'))
mv ${FULLY_RM[@]} ./exclude/


#'
#' I'm keeping all the other genes, but if you want to exclude others that
#' have at least one species with >50% removed, you could use this:
#'
# PART_RM=($(python3 << EOF
# import glob
# warns = glob.glob("*.warning")
# exclude = []
# for w in warns:
#     f = open(w, "r")
#     ll = [x.rstrip() for x in f.readlines() if "sequence removed/missing" in x]
#     f.close()
#     if len(ll) == 0:
#         continue
#     xx = [float(x.split(" ")[1].replace("%", "")) for x in ll]
#     if max(xx) > 50:
#         exclude.append(w.replace(".faa.warning", "*"))
# print(" ".join(exclude))
# EOF
# ))
# mv ${PART_RM[@]} ./exclude/


cd ..


rm -r common_genes one_filter.sh

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/

rm -r ${OUT_DIR}
