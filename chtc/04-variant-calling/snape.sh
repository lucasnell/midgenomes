#!/bin/bash

#' This script calls snape-pooled to estimate allele frequencies for
#' Pool-seq data.


#' Note that for many of the scripts inside `snape-files`, the `output`
#' argument is supposed to be a prefix only.
#' Hence, there are multiple instances below where I remove the extension
#' to the final file names for input to this argument.


# Check previous command's exit status.
# If != 0, then archive working dir and exit.
check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2" 1>&2
    cd ..
    tar -czf ERROR_${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ERROR_${OUT_DIR}.tar.gz /staging/lnell/dna/snape/
    rm -r ${OUT_DIR}
    exit $2
  fi
  echo "Checked step $1"
}


. /app/.bashrc
conda activate main-env


# Argument from submit file gives you the base of the read FASTQ file name.
# Here, I assume that FASTQ file names are of the form
# <sample name>_S<sample number>_L002_R<read number>_001.fastq.gz
# where everything not in brackets is constant.
# This is true for both my RNA and DNA sequencing reads.
# For the file `Ash-19_S5_L002_R2_001.fastq.gz`, the base of the read name
# would be "Ash-19_S5" because everything else can be inferred.

export READ_BASE=$1


#' The second input indicates the number of adult midges in this sample.
#'
export N_ADULTS=$2

# export READ_BASE=Blik-19_S6
# export N_ADULTS=40


#' ========================================================================
#' Inputs
#' ========================================================================

export IN_MP=${READ_BASE}_mpileup.txt
export GENOME=tany_scaffolds.fasta
if [ ! -f /staging/lnell/dna/mpileup/${IN_MP}.gz ]; then
    echo "/staging/lnell/dna/mpileup/${IN_MP}.gz does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 222
fi



#' ========================================================================
#' Outputs
#' ========================================================================

# Where to send everything when done:
export TARGET=/staging/lnell/dna/snape
# Final files / directories
export OUT_DIR=${READ_BASE}_snape
export SNAPE_FILE=${READ_BASE}_snape.txt.gz
export SYNC_FILE=${READ_BASE}_snape.sync.gz
export MASK_FILES_PREFIX=${READ_BASE}_snape_masked
# Intermediates
export MP_INFO_PREFIX=${IN_MP/.txt/}_info


#' ========================================================================
#' Do the things
#' ========================================================================

#' ----------------------------------------------------
#' Prep for downstream steps.

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/mpileup/${IN_MP}.gz ./ && gunzip ${IN_MP}.gz

# Files to process SNAPE output:
cp /staging/lnell/snape-files.tar.gz ./
tar -xzf snape-files.tar.gz
rm snape-files.tar.gz

# "Pickle" the reference
export GENOME_PICKLE=${GENOME/.fasta/_pickle}.ref
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz
python3 ./snape-files/PickleRef.py --ref ${GENOME} \
    --output ${GENOME_PICKLE/.ref/}
check_exit_status "PickleRef" $?

# Scaffold names in the same order as the reference.
# Used to organize output files since SNAPE analyzes by scaffold.
# Make sure scaffold names don't have spaces in them!
export SCAFF_NAMES=($(grep "^>" ${GENOME} | sed 's/>//g' | sed 's/\s.*$//'))

# Set parameters for snape:
export THETA=0.005
export D=0.01
export PRIOR_TYPE="informative"
export FOLD="unfolded"

# Set parameters for processing snape output:
export MAX_SNAPE=0.9
export MAX_COV=0.95
export MIN_COV=10




#' ----------------------------------------------------
#' Summarize some mpileup file info for use later
#' Produces the following files:
#' - ${MP_INFO_PREFIX}.sync.gz
#' - ${MP_INFO_PREFIX}.cov
#' - ${MP_INFO_PREFIX}.indel
python3 ./snape-files/Mpileup2Sync.py \
  --mpileup ${IN_MP} \
  --ref ${GENOME_PICKLE} \
  --output ${MP_INFO_PREFIX} \
  --base-quality-threshold 25 \
  --coding 1.8 \
  --minIndel 5

check_exit_status "Mpileup2Sync" $?




#' ----------------------------------------------------
#' Split mpileup file by scaffold
mkdir tmp
cd tmp

# This does the splitting:
awk '{if (last != $1) close(last); print >> $1; last = $1}' ../${IN_MP}

# Rename scaffold files and verify they exist in the reference
for scaff in *; do
    if [[ ! " ${SCAFF_NAMES[*]} " =~ " ${scaff} " ]]; then
        echo "Scaffold ${scaff} does not exist in reference!" 1>&2
        exit 1
    fi
    mv $scaff ${scaff}_mp.txt
done

# Now go back to main directory
mv *_mp.txt ../
cd ..
rm -r tmp
rm ${IN_MP}




#' ----------------------------------------------------
#' SNAPE-pooled

# I'm using `SCAFF_NAMES` so that the combined file is in the same order
# as the reference.
for scaff in ${SCAFF_NAMES[@]}; do
    if [ ! -f ${scaff}_mp.txt ]; then
        echo "mpileup output for scaffold '${scaff}' not found. Skipping..."
        continue
    fi
    ERR_FILE=${scaff}-${SNAPE_FILE/.txt.gz/.err}
    # Below, I redirect stderr to avoid printing many warnings about indels
    snape-pooled -nchr $(($N_ADULTS*2)) -theta ${THETA} -D ${D} \
        -priortype ${PRIOR_TYPE} -fold ${FOLD} < ${scaff}_mp.txt 1> \
        ${scaff}-${SNAPE_FILE/.gz/} 2> \
        ${ERR_FILE}
    check_exit_status "SNAPE-pooled (${scaff})" $?
    gzip -c ${scaff}-${SNAPE_FILE/.gz/} >> ${SNAPE_FILE}
    check_exit_status "gzip -c SNAPE output (${scaff})" $?
    rm ${scaff}_mp.txt
    # Check *.err file for anything unusual; delete if it doesn't.
    NON_WARNS=$(tail -n +5 ${ERR_FILE} | grep -v "^warning:" | wc -l)
    if [ $NON_WARNS -ge 1 ]; then
        echo -e "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" 1>&2
        echo "SNAPE call on '${scaff}' may contain errors!" 1>&2
        echo "see ${ERR_FILE}" 1>&2
        echo -e "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n" 1>&2
    else
        rm ${ERR_FILE}
    fi
done



rm *-${SNAPE_FILE/.gz/}




#' ----------------------------------------------------
#' Convert to gSYNC file
python3 ./snape-files/SNAPE2SYNC.py \
    --input ${SNAPE_FILE} \
    --ref ${GENOME_PICKLE} \
    --output ${SYNC_FILE/.sync.gz/} > \
    ${SYNC_FILE/.sync.gz/_2SYNC.out}

check_exit_status "SNAPE2SYNC" $?

# I'm capturing output above to avoid unnecessary warnings printed to stdout
# Remove unnecessary info, print file to stdout, then delete file
grep -v "^Reference is N" ${SYNC_FILE/.sync.gz/_2SYNC.out} | \
    sed '/^[[:space:]]*$/d'
rm ${SYNC_FILE/.sync.gz/_2SYNC.out}

rm ${GENOME_PICKLE}
rm ${GENOME}


#' ----------------------------------------------------
#' Create masked gSYNC file
#' Creates two files:
#' - ${MASK_FILES_PREFIX}.bed.gz
#' - ${MASK_FILES_PREFIX}_masked.sync.gz

python3 ./snape-files/MaskSYNC_snape.py \
    --sync ${SYNC_FILE} \
    --output ${MASK_FILES_PREFIX} \
    --indel ${MP_INFO_PREFIX}.indel \
    --coverage ${MP_INFO_PREFIX}.cov \
    --mincov ${MIN_COV} \
    --maxcov ${MAX_COV} \
    --maxsnape ${MAX_SNAPE} \
    --SNAPE

check_exit_status "MaskSYNC_snape" $?

# Add the following to the call above if you have a *.gff/*.gff3 of TEs:
#     --te ${TE}.gff or ${TE}.gff3 \

# This removes redundant _masked ending to this file:
mv ${MASK_FILES_PREFIX}_masked.sync.gz ${MASK_FILES_PREFIX}.sync.gz

rm ${MP_INFO_PREFIX}*

rm -r snape-files


#' ----------------------------------------------------
#' Handle output files

mv ${MASK_FILES_PREFIX}.sync.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

