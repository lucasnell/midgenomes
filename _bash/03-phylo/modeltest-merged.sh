#!/bin/bash



#'
#' Model selection using ModelTest-NG on merged partitions output from
#' ModelFinder.
#'
#'
#' Outputs:
#' - chir_modeltest_merged.tar.gz
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)


export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_modeltest_merged
export PREFIX=${OUT_DIR}
export OUT_PART=${PREFIX}.partition

mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz

# --------
# merged partitions:
export ALIGNS_PART=chir_merge.partition
cp ${TARGET}/${ALIGNS_PART} ./

# Replace all models with 'AA'
sed -i 's/^[^,]*,/AA,/' ${ALIGNS_PART}


#'
#' -t ml
#' Optimize topology and branch lengths for each model
#'
#' -d aa
#' amino acid
#'
#' -m DAYHOFF,LG,DCMUT,JTT,WAG,BLOSUM62,JTT-DCMUT,MTREV,MTART,MTZOA
#' sets the candidate model matrices separated by commas
#' I chose all models that are for nuclear or mitochondrial sources,
#'  where MT sources are filtered for those relevant to dipterans.
#' I also removed GTR because it takes too long to fit in RAxML-NG,
#' and I removed VT and PMB because they aren't available in MCMCTree.
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
#' -h g
#' sets the candidate models rate heterogeneity to include
#' discrete Gamma rate categories (+G)
#'

modeltest-ng -i ${CONCAT_ALIGNS} \
    -t ml \
    -d aa \
    -m DAYHOFF,LG,DCMUT,JTT,WAG,BLOSUM62,JTT-DCMUT,MTREV,MTART,MTZOA \
    -o ${PREFIX} \
    -q ${ALIGNS_PART} \
    -p ${THREADS} \
    -r 837645180 \
    -h g

# These are stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}


cp ${PREFIX}.part.aic ../${OUT_PART}

cd ..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_PART} ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

