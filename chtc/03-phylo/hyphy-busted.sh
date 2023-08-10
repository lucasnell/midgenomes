#!/bin/bash


#'
#' Selection tests using HyPhy BUSTED method.
#'


. /app/.bashrc

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_hyphy_busted
mkdir ${OUT_DIR}
cd ${OUT_DIR}

#' Ultrametric species tree without Cmarin (bc of poor matching in OrthoFinder)
export SPECIES_TREE=chir_mcmctree_noCmarin.nwk

cp ${TARGET}/phylo/${SPECIES_TREE} ./
check_exit_status "moving species tree" $?


export ALIGNS_DIR=hyphy_aligns
export ALIGNS_PARENT_DIR=chir_hyphy_inputs
tar -xf ${TARGET}/${ALIGNS_PARENT_DIR}.tar.gz -C ./

cd ${ALIGNS_PARENT_DIR}
mv ${ALIGNS_DIR} ../
cd ..
rm -r ${ALIGNS_PARENT_DIR}

cd ${ALIGNS_DIR}
export HOG_NAMES=($(find . -type d -name "N0*" | sed 's/\.\///g'))
cd ..


# Script to run `hyphy busted` on one HOG:
cat << EOF > one_hyphy.sh
#!/bin/bash
hog=\$1

align=./${ALIGNS_DIR}/\${hog}/\${hog}.msa
hyphy CPU=1 busted --alignment \${align} \\
    --tree ${SPECIES_TREE} \\
    --output \${hog}.json \\
    --save-fit \${hog}.fit \\
    1> \${hog}.stdout \\
    2> \${hog}.stderr
status=\$?
exit \$status
EOF
chmod +x one_hyphy.sh



python3 << EOF
import subprocess as sp
import multiprocessing as mp
import sys

def work(hog):
    """Defines the work unit on an input file"""
    cmd = "./one_hyphy.sh " + hog
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    tasks = "${HOG_NAMES[@]}".split(" ")
    n_tasks = len(tasks)
    with mp.Pool(processes=${THREADS}) as pool:
        return_codes = pool.map(work, tasks)
        return_codes = [x.returncode for x in return_codes]
        print("-----------------\nreturn codes:\n-----------------")
        for i in range(n_tasks):
            print(tasks[i] + " = " + str(return_codes[i]))
    sys.exit(0)

EOF



rm -r ${ALIGNS_DIR}

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}
