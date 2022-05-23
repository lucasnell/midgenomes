#!/bin/bash

# Use purge_dups on NECAT and SMARTdenovo assemblies

. /app/.bashrc
conda activate assembly-env
source /staging/lnell/helpers.sh

# Where to send files:
export TARGET=/staging/lnell/assemblies

export ASSEMBLIES=(contigs_necat_ont-polish_nextpolish.fasta \
                   contigs_smart_ont-polish_nextpolish.fasta)
export L_CUTS=(30 40)
export M_CUTS=(150 155)
export U_CUTS=(390 400)
export a_ARGS=(70 45)


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

# Make temporary working directory to delete later:
export WD=work_dir
mkdir ${WD}
cd ${WD}


# Initial alignments

export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/${LONGREADS} ./

for ASS in ${ASSEMBLIES[@]}; do
    cp /staging/lnell/assemblies/${ASS}.gz ./ && gunzip ${ASS}.gz
    ALIGN=${ASS/.fasta/_mm2-align.paf.gz}
    minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 ${ASS} ${LONGREADS} | \
        gzip -c - > ${ALIGN}
done

rm ${LONGREADS}


for i in 0 1; do

    A=${ASSEMBLIES[i]/.fasta/}
    OUT_DIR=${A}_purgedups
    mkdir ${OUT_DIR}
    cd ${OUT_DIR}
    ALIGN=../${A}_mm2-align.paf.gz
    OUT_FASTA=${OUT_DIR}.fasta
    #  (produces PB.base.cov and PB.stat files)
    pbcstat ${ALIGN}
    # calculate cutoffs, setting limits manually
    calcuts -l ${L_CUTS[i]} -m ${M_CUTS[i]} -u ${U_CUTS[i]} PB.stat \
        > cutoffs \
        2> calcuts.log

    # Split an assembly and do a self-self alignment:
    split_fa ../${A}.fasta > ${A}.fasta.split
    minimap2 -x asm5 -DP ${A}.fasta.split ${A}.fasta.split \
        -t $((THREADS - 2)) -K 1G -2 | \
        gzip -c - > ${A}.fasta.split.self.paf.gz

    # Purge haplotigs and overlaps:
    purge_dups -2 -T cutoffs -c PB.base.cov \
        -a ${a_ARGS[i]} \
        ${A}.fasta.split.self.paf.gz > \
        dups.bed 2> purge_dups.log

    # Get purged primary and haplotig sequences from draft assembly:
    get_seqs -e dups.bed ../${A}.fasta

    mv purged.fa ${OUT_FASTA}

    summ-scaffs.py ${OUT_FASTA} | \
        tee scaff_summary.out
    run_busco ${OUT_FASTA} ${THREADS}
    rm -r busco busco_downloads
    busco_seq_summary_csv scaff_summary.out busco.out ${OUT_FASTA/.fasta/}

    gzip ${OUT_FASTA}
    mv ${OUT_FASTA}.gz ${TARGET}/
    mv ${ALIGN} ./
    cd ..
    tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ${OUT_DIR}.tar.gz ${TARGET}/
    rm -r ${OUT_DIR}

done

rm *.fasta

cd ..
rm -r ${WD}

