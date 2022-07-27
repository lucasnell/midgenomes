#!/bin/bash


#'
#' Run 100 bootstrap replicates using RAxML-NG.
#' It requires as input an index for which seed to use (0-9).
#' This allows bootstrapping to be split among multiple smaller jobs to be
#' combined later.
#'


if (( $# != 1 )); then
    echo -e "ERROR: raxml-boot.sh requires exactly 1 argument" 1>&2
    exit 1
fi
export SEEDS_IND=$1
if ! [[ $SEEDS_IND =~ ^[0-9]+$ ]] || (( SEEDS_IND < 0 )) || (( SEEDS_IND > 9 )); then
    echo -n "ERROR: arg to raxml-boot.sh should be an integer 0-9. " 1>&2
    echo "Yours is '${SEEDS_IND}'." 1>&2
    exit 2
fi
#' Array of seeds (previously generated using `sample.int(2^31-1, 10)` in R):
SEEDS_ARR=(1686233904  838153601 1341023619 1137233470 1011779318 \
           1309257016 1615251508 1755012526 1662390602 1658268201)

export THIS_SEED=${SEEDS_ARR[$SEEDS_IND]}


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_boot_${SEEDS_IND}
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}



export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
check_exit_status "move, extract alignments" $?



raxml-ng --bootstrap --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} --threads ${THREADS} \
    --bs-trees 10 \
    --data-type AA \
    --model LG+I+G \
    --seed ${THIS_SEED} \
    1> >(tee -a ${PREFIX}.stdout)

##    --bs-trees 100 \
##    --outgroup Anopheles_stephensi \

rm ${CONCAT_ALIGNS}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



