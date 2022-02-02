#!/bin/bash

# Explore options for purge_dups

. /app/.bashrc
conda activate main-env


export THREADS=24

export REF=contigs_shasta_pepper-hap1.fasta
cp /staging/lnell/${REF}.gz ./ && gunzip ${REF}.gz


# export THREADS=16
# export LONGREADS=basecalls_guppy-5.0.11.fastq
# cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
# conda activate main-env
# minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 \
#     ${REF} ${LONGREADS} | \
#     gzip -c - > ont_align_minimap2.paf.gz
#
# mv ont_align_minimap2.paf.gz /staging/lnell/
# rm ${REF} ${LONGREADS}

export ALIGN=ont_align_minimap2.paf.gz
cp /staging/lnell/${ALIGN} ./

export OUT_DIR=epd_test

mkdir ${OUT_DIR}
mv ${REF} ${ALIGN} ./${OUT_DIR}/
cd ${OUT_DIR}


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
    -l 5K \
    ${REF}.split.self.paf.gz > \
    dups.bed 2> purge_dups.log

    # -E 12K \

    # -G 100K \
    # -M 40K \
    # -m 1000 \
    # -b 500 \
    # -f 0.9 \
    # -a 90 \

# -f  INT   minimum fraction of haploid/diploid/bad/repetitive bases in a sequence [.8]
# -a  INT   minimum alignment score [70]
# -b  INT   minimum max match score [200]
# -2  BOOL  2 rounds chaining [FALSE]
# -m  INT   minimum matching bases for chaining [500]
# -M  INT   maximum gap size for chaining [20K]
# -G  INT   maximum gap size for 2nd round chaining [50K]
# -l  INT   minimum chaining score for a match [10K]
# -E  INT   maximum extension for contig ends [15K]


# Get purged primary and haplotig sequences from draft assembly:
get_seqs -e dups.bed ${REF}

# export OUT_FASTA=purged_${REF}
# mv purged.fa ${OUT_FASTA}
export OUT_FASTA=purged.fa

summ-scaffs.py ${OUT_FASTA}

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS}
conda deactivate

rm -r busco ${OUT_FASTA} dups.bed purge_dups.log hap.fa

chmod +x pdi.sh


sed "s/__ARGS__/-a\ 93/g" ../pd.sh > pdi.sh
./pdi.sh






for i in 92 94 95
do
    echo -e "\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
    echo -e ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
    echo -e " a = ${i}\n"
    echo -e ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
    echo -e ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n"
    sed "s/__ARGS__/-a\ ${i}/g" ../pd.sh > pdi.sh
    ./pdi.sh
done



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}

# # If you want to use anything in the scripts folder,
# # it's not in the conda package.
# wget https://github.com/dfguan/purge_dups/archive/refs/tags/v1.2.5.tar.gz
# tar -xzf v1.2.5.tar.gz
# rm v1.2.5.tar.gz
# cd purge_dups-1.2.5
# mv scripts ../
# cd ..
# rm -r purge_dups-1.2.5
#
#
# scripts/hist_plot.py




# Another file `pd.sh` that does just the purge_dups and onward steps
# to look at different args to purge_dups


#!/bin/bash


. /app/.bashrc
conda activate main-env


export REF=contigs_shasta_pepper-hap1.fasta
export THREADS=24
export OUT_FASTA=purged.fa

# Purge haplotigs and overlaps:
purge_dups -2 -T cutoffs -c PB.base.cov \
    __ARGS__ \
    ${REF}.split.self.paf.gz > \
    dups.bed 2> purge_dups.log

# Get purged primary and haplotig sequences from draft assembly:
get_seqs -e dups.bed ${REF}

summ-scaffs.py ${OUT_FASTA}

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS}
conda deactivate

rm -r busco ${OUT_FASTA} dups.bed purge_dups.log hap.fa
