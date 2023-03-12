#!/bin/bash

#' Use purge_dups on NECAT and SMARTdenovo assemblies.
#' This does not come with a *.sub file because it was run interactively.

. /app/.bashrc
conda activate assembly-env

# Where to send files:
export TARGET=/staging/lnell/assemblies


export ASSEMBLIES=(contigs_necat_ont-polish.fasta \
                   contigs_smart_ont-polish.fasta \
                   contigs_flye_ont-polish.fasta)
export L_CUTS=(10 30 10)
export M_CUTS=(150 160 150)
export U_CUTS=(390 400 400)
export a_ARGS=(70 45 60)


export THREADS=$(count_threads)

# Make temporary working directory to delete later:
export WD=work_dir
mkdir ${WD}
cd ${WD}


# Initial alignments

export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/ont/${LONGREADS} ./

export N_ASS=${#ASSEMBLIES[@]}

for ASS in ${ASSEMBLIES[@]}; do
    cp /staging/lnell/assemblies/${ASS}.gz ./ && gunzip ${ASS}.gz
    ALIGN=${ASS/.fasta/_mm2-align.paf.gz}
    minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 ${ASS} ${LONGREADS} | \
        gzip -c - > ${ALIGN}
done

rm ${LONGREADS}


for ((i=0; i < N_ASS; i+=1)); do
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
        tee contigs_summary.out
    run_busco ${OUT_FASTA} ${THREADS}
    rm -r busco busco_downloads
    pretty-csv.py -s contigs_summary.out -b busco.out ${OUT_FASTA/.fasta/}

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

