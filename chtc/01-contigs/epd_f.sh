#!/bin/bash

# Explore options for purge_dups -f argument
# usage: ./epd.sh [VALUE FOR `purge_dups -f` ARG]


. /app/.bashrc
conda activate main-env

export REF=contigs_shasta_pepper-hap1.fasta
export THREADS=12

export F_ARG_VAL=$1

export OUT_PREFIX=purgedups_f${F_ARG_VAL}
export OUT_DIR=${OUT_PREFIX}
export OUT_FASTA=${OUT_PREFIX}.fasta
export OUT_CSV=${OUT_PREFIX}.csv
export OUT_PNG=${OUT_PREFIX}.png

# Make temporary working directory to delete later:
export WD=work_dir
mkdir ${WD}
cd ${WD}

cp /staging/lnell/${REF}.gz ./ && gunzip ${REF}.gz

# Initial alignment run separately using the following code:
# export LONGREADS=basecalls_guppy-5.0.11.fastq
# export OUT_PAF=ont_align_mm2_pepper.paf.gz
# cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
# conda activate main-env
# minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 \
#     ${REF} ${LONGREADS} | \
#     gzip -c - > ${OUT_PAF}
#
# mv ${OUT_PAF} /staging/lnell/
# rm ${REF} ${LONGREADS}

export ALIGN=ont_align_mm2_pepper.paf.gz
cp /staging/lnell/${ALIGN} ./

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
# good one.
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


summ-scaffs.py ${OUT_FASTA} > \
    scaff_summary.out

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS} > \
    busco.out
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
export LONGREADS=basecalls_guppy-5.0.11.fastq
export NEW_ALIGN=ont_align_mm2_${OUT_PREFIX}.paf.gz
cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
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
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}

cd ..
rm -r ${WD}

