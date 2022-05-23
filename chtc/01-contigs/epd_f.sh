#!/bin/bash

#'
#' Explore options for purge_dups -a and -f arguments
#' usage: ./epd.sh [REFERENCE ASSEMBLY] \
#'                 [ALIGNMENT] \
#'                 [VALUE FOR `purge_dups -a` ARG] \
#'                 [VALUE FOR `purge_dups -f` ARG]
#'



. /app/.bashrc
conda activate assembly-env

export REF=$1
export ALIGN=$2
export A_ARG_VAL=$3
export F_ARG_VAL=$4

export LONGREADS=basecalls_guppy-5.0.11.fastq
export THREADS=12

# Check for correct suffixes:
if [[ ${REF} != *.fasta ]]; then
    echo -e "\n\nERROR: ${REF} doesn't end in '.fasta'" 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi
if [[ ${ALIGN} != *.paf.gz ]]; then
    echo -e "\n\nERROR: ${ALIGN} doesn't end in '.paf.gz'" 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi

export TARGET=/staging/lnell/epd
export OUT_PREFIX=purgedups_a-${A_ARG_VAL}_f-${F_ARG_VAL}
export OUT_DIR=${OUT_PREFIX}
export OUT_FASTA=${OUT_PREFIX}.fasta
export OUT_CSV=${OUT_PREFIX}.csv
export OUT_PNG=${OUT_PREFIX}.png

# Make temporary working directory to delete later:
export WD=work_dir
mkdir ${WD}
cd ${WD}

# Check previous command's exit status.
# If != 0, then archive working dir and exit.
check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2" 1>&2
    cd ..
    mv ${WD} ERROR_${OUT_PREFIX}
    tar -czf ERROR_${OUT_PREFIX}.tar.gz ERROR_${OUT_PREFIX}
    mv ERROR_${OUT_PREFIX}.tar.gz ${TARGET}/
    rm -r ERROR_${OUT_PREFIX}
    exit $2
  fi
  echo "Checked step $1"
}

cp /staging/lnell/${REF}.gz ./ && gunzip ${REF}.gz
cp /staging/lnell/${ALIGN} ./
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz




# Initial alignment run separately using the following code:
# export OUT_PAF=ont_align_mm2_pepper.paf.gz
# cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
# conda activate main-env
# minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 \
#     ${REF} ${LONGREADS} | \
#     gzip -c - > ${OUT_PAF}
#
# mv ${OUT_PAF} /staging/lnell/
# rm ${REF} ${LONGREADS}





# ======================================================================
# ======================================================================
# Initial purging
# ======================================================================
# ======================================================================

#  (produces PB.base.cov and PB.stat files)
pbcstat ${ALIGN}
# calculate cutoffs, setting upper limit manually
# purge_dups creators say highly heterozygous genomes can have too-high
# upper limit, and based on histogram from hist_plot.py, 380 seems like a
# good one here.
calcuts -l 5 -m 181 -u 380 PB.stat > cutoffs 2> calcults.log
# originally:
# 5	91	151	181	302	543
# now:
# 5	180	180	181	181	380

# Split an assembly and do a self-self alignment:
split_fa ${REF} > ${REF}.split
minimap2 -x asm5 -DP ${REF}.split ${REF}.split \
    -t $((THREADS - 2)) -K 1G -2 | \
    gzip -c - > ${REF}.split.self.paf.gz


# Purge haplotigs and overlaps:
purge_dups -2 -T cutoffs -c PB.base.cov \
    -f ${F_ARG_VAL} \
    ${REF}.split.self.paf.gz > \
    dups.bed 2> purge_dups.log

# Get purged primary and haplotig sequences from draft assembly:
get_seqs -e dups.bed ${REF}

mv purged.fa ${OUT_FASTA}

# ======================================================================
# ======================================================================
# Summarize initial purged assembly
# ======================================================================
# ======================================================================


summ-scaffs.py ${OUT_FASTA} | \
    tee scaff_summary.out

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS} | \
    tee busco.out
conda deactivate

# output args identifier:
echo -n "${OUT_PREFIX}", >> ${OUT_CSV}
# output from summ-scaffs.py:
echo -n $(grep "size" scaff_summary.out | sed 's/.* //'), >> ${OUT_CSV}
echo -n $(grep "scaffolds$" scaff_summary.out | sed -r 's/\ .+//'), >> ${OUT_CSV}
echo -n $(grep "N50" scaff_summary.out | sed 's/.* //'), >> ${OUT_CSV}
echo -n $(grep "min" scaff_summary.out | sed 's/.* //'), >> ${OUT_CSV}
echo -n $(grep "max" scaff_summary.out | sed 's/.* //'), >> ${OUT_CSV}
echo -n $(grep "total\ N" scaff_summary.out | sed 's/.* //'), >> ${OUT_CSV}
# output from BUSCO:
busco_var() {
    grep "|" busco.out | grep "${BVS}" | tr -d "|" | xargs | sed 's/ .*//'
}
echo -n $(BVS="(C)" && busco_var), >> ${OUT_CSV}
echo -n $(BVS="(S)" && busco_var), >> ${OUT_CSV}
echo -n $(BVS="(D)" && busco_var), >> ${OUT_CSV}
echo -n $(BVS="(F)" && busco_var), >> ${OUT_CSV}
echo -n $(BVS="(M)" && busco_var), >> ${OUT_CSV}
echo $(BVS="Total BUSCO" && busco_var) >> ${OUT_CSV}

cat ${OUT_CSV}


# ======================================================================
# ======================================================================
# Make read depth histogram plot of purged assembly
# ======================================================================
# ======================================================================

# Align ONT reads to purged genome
export NEW_ALIGN=ont_align_mm2_${OUT_PREFIX}.paf.gz
minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 \
    ${OUT_FASTA} ${LONGREADS} | \
    gzip -c - > ${NEW_ALIGN}

# Make new PB.stat file
pbcstat ${NEW_ALIGN}

# If you want to use anything in purge_dups's scripts folder,
# it's not in the conda package.
wget https://github.com/dfguan/purge_dups/archive/refs/tags/v1.2.5.tar.gz
tar -xzf v1.2.5.tar.gz
rm v1.2.5.tar.gz

# Make histogram plot
purge_dups-1.2.5/scripts/hist_plot.py PB.stat ${OUT_PNG}


mkdir ${OUT_DIR}
mv ${OUT_FASTA} ${OUT_CSV} ${OUT_PNG} ./${OUT_DIR}/
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}

cd ..
rm -r ${WD}

