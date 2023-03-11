#!/bin/bash


#'
#' Pre-alignment quality filter with prequal
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

conda activate phylo-env

mkdir working
cd working

cp /staging/lnell/phylo/odb/shared_odb.tar.gz ./
tar -xzf shared_odb.tar.gz
rm shared_odb.tar.gz


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


#' Only show progess bar below if in an interactive job:
export inter_job="False"
if [[ $- == *i* ]]; then inter_job="True"; fi


python3 << EOF
import glob
import subprocess as sp
import multiprocessing as mp
import sys

def work(fa_file):
    """Defines the work unit on an input file"""
    cmd = "./one_filter.sh " + fa_file
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    tasks = glob.glob("shared_odb/*.faa")
    n_tasks = len(tasks)
    show_progress = $inter_job
    with mp.Pool(processes=${THREADS}) as pool:
        for i, _ in enumerate(pool.imap_unordered(work, tasks), 1):
            if show_progress:
                sys.stderr.write('\rdone {0:%}'.format(i/n_tasks))
    if show_progress:
        sys.stderr.write('\n')
    sys.exit(0)

EOF




#'
#' Exclude any genes where one or more species had sequences fully removed
#'
cd ${OUT_DIR}
FULLY_RM=($(egrep -lir --include=*.warning "Fully removed" . \
    | sed 's/\.faa\.warning/\*/g'))
if (( ${#FULLY_RM[@]} > 0 )); then
    mkdir exclude
    mv ${FULLY_RM[@]} ./exclude/
fi


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


rm -r shared_odb one_filter.sh

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/

cd ..
rm -r working
