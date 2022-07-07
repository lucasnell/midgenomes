#!/bin/bash


#'
#' Sequence alignment using MAFFT.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

conda activate phylo-env

cp /staging/lnell/phylo/prequal_filt.tar.gz ./
tar -xzf prequal_filt.tar.gz
rm prequal_filt.tar.gz


export OUT_DIR=mafft_aligns
mkdir ${OUT_DIR}




#' Script to run mafft to do alignment on a single gene.
#' Usage:
#'      ./one_align.sh [*.faa FILE]
cat << EOF > one_align.sh
#!/bin/bash
IN=\$1
NAME_BASE=\$(basename \$1 | sed 's/\..*//g')
OUT=./${OUT_DIR}/\${NAME_BASE}.faa
LOG=./${OUT_DIR}/\${NAME_BASE}.log
linsi --thread 1 \${IN} > \${OUT} 2> \${LOG}
EOF
chmod +x one_align.sh


python3 << EOF
import glob
import subprocess as sp
import multiprocessing as mp
import sys

def work(fa_file):
    """Defines the work unit on an input file"""
    cmd = "./one_align.sh " + fa_file
    ret = sp.run(cmd, shell = True)
    return ret

tasks = glob.glob("prequal_filt/*.faa.filtered")
n_tasks = len(tasks)

with mp.Pool(processes=${THREADS}) as pool:
    for i, _ in enumerate(pool.imap_unordered(work, tasks), 1):
        sys.stderr.write('\rdone {0:%}'.format(i/n_tasks))
sys.stderr.write('\n')

sys.exit(0)

EOF


rm -r prequal_filt one_align.sh

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/


rm -r ${OUT_DIR}
