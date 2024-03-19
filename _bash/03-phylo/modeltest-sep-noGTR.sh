#!/bin/bash


#'
#' Model selection using ModelTest-NG.
#'
#'
#' Outputs:
#' - chir_modeltest_noGTR_${MT_INDEX}.tar.gz
#'

if (( $# != 1 )); then
    echo "ERROR: Exactly one argument required for modeltest.sh" 1>&2
    exit 1
fi

export MT_INDEX=$1

if ! [[ $MT_INDEX =~ ^[0-9]+$ ]] || (( MT_INDEX < 0 )) ||
        (( MT_INDEX > 9 )); then
    echo -n "ERROR: First argument should be an integer 0-9. " 1>&2
    echo "Yours is '${MT_INDEX}'." 1>&2
    exit 1
fi

RNG_SEEDS=(800899873 269997286 362470401 61320647 996066934 \
           1085297174 221024413 1350768726 23079768 1075989655)
export RNG_SEED=${RNG_SEEDS[$MT_INDEX]}
unset -v RNG_SEEDS


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)


export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_modeltest_noGTR_${MT_INDEX}
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz

# --------
# unmerged partitions:
export ALIGNS_PART=trim_aligns.partition
cp ${TARGET}/${ALIGNS_PART}.gz ./ && gunzip ${ALIGNS_PART}.gz


# python code to split partition into 10 separate files
# (we don't need to do this to the alignment):
python3 << EOF
import sys
with open("${ALIGNS_PART}", "r") as file:
    parts = [x for x in file.readlines()]
starts = [round(len(parts) * x / 10) for x in range(10)]
ends = [x-1 for x in starts[1:]]
ends.append(len(parts)-1)
i=${MT_INDEX}
fn = "${ALIGNS_PART%.*}_" + str(i) + ".partition"
n_lines = ends[i] - starts[i] + 1
with open(fn, "w") as f:
    for j in range(n_lines):
        f.write(parts[j+starts[i]])
sys.exit(0)
EOF

rm ${ALIGNS_PART}

#'
#' -t ml
#' Optimize topology and branch lengths for each model
#'
#' -d aa
#' amino acid
#'
#' -m DAYHOFF,LG,DCMUT,JTT,WAG,VT,BLOSUM62,PMB,JTT-DCMUT,MTREV,MTART,MTZOA
#' sets the candidate model matrices separated by commas
#' I chose all models that are for nuclear or mitochondrial sources,
#'  where MT sources are filtered for those relevant to dipterans.
#' I also removed GTR because it takes too long to fit in RAxML-NG.
#'
#' -o modeltest_${MT_INDEX}.log
#' pipes the output into a file
#'
#' -q ${ALIGNS_PART}
#' sets partitions file
#'
#' -p ${THREADS}
#' number of processes to use (shared memory)
#'
#' -r INTEGER
#' sets the seed for the random number generator
#'
#' -h f
#' sets the candidate models rate heterogeneity to include both
#' proportion of invariant sites (+I) and discrete Gamma rate categories (+G)
#'

modeltest-ng -i ${CONCAT_ALIGNS} \
    -t ml \
    -d aa \
    -m DAYHOFF,LG,DCMUT,JTT,WAG,VT,BLOSUM62,PMB,JTT-DCMUT,MTREV,MTART,MTZOA \
    -o modeltest_${MT_INDEX}.log \
    -q ${ALIGNS_PART%.*}_${MT_INDEX}.partition \
    -p ${THREADS} \
    -r ${RNG_SEED} \
    -h f

# This is stored elsewhere:
rm ${CONCAT_ALIGNS}


cd ..


tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

