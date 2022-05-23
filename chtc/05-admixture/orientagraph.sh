#!/bin/bash


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

export TARGET=/staging/lnell/admixture

export N_MIG=$1

export OUT_PREFIX=oag_m${N_MIG}
export OUT_DIR=${OUT_PREFIX}
mkdir ${OUT_DIR}
cd ${OUT_DIR}

. /app/.bashrc
source /staging/lnell/helpers.sh
conda activate orientagraph-env



SNP=space_sync_masked_biallelic.snp.gz
cp /staging/lnell/admixture/${SNP} ./
export SNP=$(pwd)/${SNP}

cp /staging/lnell/admixture/oag_seeds.txt.gz ./
SEEDS=($(gunzip -c oag_seeds.txt.gz | cut -f $(( N_MIG + 1 )) ))
export SEEDS="${SEEDS[@]:0:${THREADS}}"
rm oag_seeds.txt.gz

# The seeds file above was generated with this R code:
# set.seed(834778645)
# matrix(sample.int(2^31-1, 11 * 100), 100, 11) %>%
#     as.data.frame() %>%
#     set_names(letters[1:11]) %>%
#     as_tibble() %>%
#     write_tsv("~/_data/snape/admixture/oag_seeds.txt.gz", col_names = FALSE)



cat << EOF > one_oag_run.sh
#!/bin/bash
export SEED=\${1}
export OUT_FILE_PREFIX=m${N_MIG}_s\${SEED}

orientagraph -i ${SNP} \\
    -seed \${SEED} \\
    -m ${N_MIG} \\
    -k 100 \\
    -global \\
    -allmigs -mlno \\
    -o \${OUT_FILE_PREFIX} \\
    1> \${OUT_FILE_PREFIX}.out \\
    2> \${OUT_FILE_PREFIX}.err

exit 0
EOF
chmod +x one_oag_run.sh


cat << EOF > mp_oag.py
#!/usr/bin/env python3
import subprocess as sp
import multiprocessing as mp
import sys
import argparse

def work(rng_seed):
    """Defines the work unit on one rep"""
    ret = sp.call(["./one_oag_run.sh", rng_seed])
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run npstat in parallel")
    parser.add_argument("inputs", nargs = "+", help="seeds")
    args = parser.parse_args()
    tasks = args.inputs
    with mp.Pool(processes=${THREADS}) as pool:
        exit_codes = pool.map(work, tasks)
    print("exit codes : {}\n".format(exit_codes))
    sys.exit(0)

EOF
chmod +x mp_oag.py


./mp_oag.py ${SEEDS}


for s in ${SEEDS}; do
    lf=m${N_MIG}_s${s}.llik
    LLIK=$(tail -n 1 ${lf} \
           | sed 's/.*ln(likelihood)\ //g' \
           | sed 's/\ with.*//g')
    echo "LogLik for seed ${s} = ${LLIK}"
done


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}
