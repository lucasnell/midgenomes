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
SEEDS_ARR=(1478865727 158688736 548184584 1400912585 357834478 667771158
           145164412 1949248987 532694851 634149308)

export THIS_SEED=${SEEDS_ARR[$SEEDS_IND]}


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml_boot_${SEEDS_IND}
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}



export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz
check_exit_status "move, extract alignments" $?
export ALIGNS_PART=chir_modeltest.partition
cp ${TARGET}/${ALIGNS_PART} ./
check_exit_status "move partitions" $?



#' Note that `--extra thread-pin` was added to avoid the following error:
#' "ERROR: CPU core oversubscription detected! RAxML-NG will terminate now
#'  to avoid wasting resources."
raxml-ng --bootstrap --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
    --bs-trees 100 \
    --threads ${THREADS} \
    --extra thread-pin \
    --outgroup Mdomes \
    --data-type aa \
    --model ${ALIGNS_PART} \
    --brlen scaled \
    --seed ${THIS_SEED} \
    1> >(tee -a ${PREFIX}.stdout)
status=$?


#' Sometimes we still get the 'CPU core oversubscription' error by chance,
#' so if that happens, we'll try it again with `--extra thread-nopin` to see
#' if that fixes the problem.
if [ "$status" != "0" ]; then
    if grep -Fq 'ERROR: CPU core oversubscription detected!' ${PREFIX}.raxml.log; then
        mv ${CONCAT_ALIGNS} ${ALIGNS_PART} ../ \
            && rm * \
            && mv ../${CONCAT_ALIGNS} ../${ALIGNS_PART} ./
        raxml-ng --bootstrap --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} \
            --bs-trees 100 \
            --threads ${THREADS} \
            --extra thread-nopin \
            --outgroup Mdomes \
            --data-type aa \
            --model ${ALIGNS_PART} \
            --brlen scaled \
            --seed ${THIS_SEED} \
            1> >(tee -a ${PREFIX}.stdout)
        check_exit_status "raxml-ng thread-nopin" $?
    else
        check_exit_status "raxml-ng" $status
    fi
fi



# These are stored elsewhere:
rm ${CONCAT_ALIGNS} ${ALIGNS_PART}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



