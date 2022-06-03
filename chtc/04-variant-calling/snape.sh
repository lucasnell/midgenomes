#!/bin/bash

#' This script calls snape-pooled to estimate allele frequencies for
#' Pool-seq data.


#' Note that for many of the scripts inside `snape-files`, the `output`
#' argument is supposed to be a prefix only.
#' Hence, there are multiple instances below where I remove the extension
#' to the final file names for input to this argument.


. /app/.bashrc
conda activate main-env
source /staging/lnell/helpers.sh


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



#' ========================================================================
#' Inputs
#' ========================================================================

export IN_MP=${READ_BASE}_mpileup.txt
export GENOME=tany_scaffolds.fasta
export REPEATS=tany_repeat_anno.gff3



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
export PART_MASK_FILES_PREFIX=${READ_BASE}_snape_part_masked
# Intermediates
export MP_INFO_PREFIX=${IN_MP/.txt/}_info




#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/mpileup/${IN_MP}.gz ./ && gunzip ${IN_MP}.gz
check_exit_status "cp ${IN_MP}.gz" $?
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp ${GENOME}.gz" $?
cp /staging/lnell/${REPEATS}.gz ./ && gunzip ${REPEATS}.gz
check_exit_status "cp ${REPEATS}.gz" $?

# Files to process SNAPE output:
cp /staging/lnell/snape-files.tar.gz ./
check_exit_status "cp snape-files.tar.gz" $?
tar -xzf snape-files.tar.gz
rm snape-files.tar.gz

# "Pickle" the reference
export GENOME_PICKLE=${GENOME/.fasta/_pickle}.ref
python3 ./snape-files/PickleRef.py --ref ${GENOME} \
    --output ${GENOME_PICKLE/.ref/}
check_exit_status "PickleRef" $?

# Scaffold names in the same order as the reference.
# Used to organize output files since SNAPE analyzes by scaffold.
# Make sure scaffold names don't have spaces in them!
export SCAFF_NAMES=($(grep "^>" ${GENOME} | sed 's/>//g' | sed 's/\s.*$//'))

# Set parameters for snape:
export THETA=0.008
export D=0.01
export PRIOR_TYPE="informative"
export FOLD="unfolded"

# Set parameters for processing snape output:
export MAX_SNAPE=0.9
export MAX_COV=0.95
export MIN_COV=10




#' ========================================================================
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




#' ========================================================================
#' Split mpileup file by scaffold
mkdir tmp
cd tmp

# This does the splitting:
awk '{if (last != $1) close(last); print >> $1; last = $1}' ../${IN_MP}

# Rename scaffold files and verify they exist in the reference
for scaff in *; do
    if [[ ! " ${SCAFF_NAMES[*]} " =~ " ${scaff} " ]]; then
        echo "Scaffold ${scaff} does not exist in reference! " 1>&2
        exit 1
    fi
    mv $scaff ${scaff}_mp.txt
done

# Now go back to main directory
mv *_mp.txt ../
cd ..
rm -r tmp
rm ${IN_MP}




#' ========================================================================
#' SNAPE-pooled

# I'm using `SCAFF_NAMES` so that the combined file is in the same order
# as the reference.
for scaff in ${SCAFF_NAMES[@]}; do
    # If this scaffold isn't present, we add a filler line to the snape
    # output file to make sure that all of our eventual gSYNC files have the
    # same number of rows.
    if [ ! -f ${scaff}_mp.txt ]; then
        echo "mpileup output for scaffold '${scaff}' not found. " \
            "Adding filler and skipping..."
        NT=$(grep -A1 "^>${scaff}$" ${GENOME} | tail -n 1 | head -c 1)
        echo -e "${scaff}\t1\t${NT}\t1\t0\t1\t1\t${NT}\t0.0\t0.0\t0.0" | \
            gzip \
            >> ${SNAPE_FILE}
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
        echo "SNAPE call on '${scaff}' may contain errors! " 1>&2
        echo "see ${ERR_FILE}" 1>&2
        echo -e "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n" 1>&2
    else
        rm ${ERR_FILE}
    fi
done





rm *-${SNAPE_FILE/.gz/}




#' ========================================================================
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


#' ========================================================================
#' Create masked gSYNC file
#' Creates two files:
#' - ${MASK_FILES_PREFIX}.bed.gz
#' - ${MASK_FILES_PREFIX}_masked.sync.gz

python3 ./snape-files/MaskSYNC_snape.py \
    --sync ${SYNC_FILE} \
    --output ${MASK_FILES_PREFIX} \
    --indel ${MP_INFO_PREFIX}.indel \
    --coverage ${MP_INFO_PREFIX}.cov \
    --te ${REPEATS} \
    --mincov ${MIN_COV} \
    --maxcov ${MAX_COV} \
    --maxsnape ${MAX_SNAPE} \
    --SNAPE

check_exit_status "MaskSYNC_snape" $?


# This removes redundant _masked ending to this file:
mv ${MASK_FILES_PREFIX}_masked.sync.gz ${MASK_FILES_PREFIX}.sync.gz

# Create a snp file for use in npstat:
gunzip -c ${MASK_FILES_PREFIX}.sync.gz \
    | grep -v ".:.:.:.:.:." \
    | cut -f 1,2 \
    | gzip \
    > ${MASK_FILES_PREFIX}.snp.gz



#' ========================================================================
#' Create partially masked gSYNC file for npstat
#' Differs from above bc it doesn't include filtering based on P(polymorphic)
#' Creates two files:
#' - ${PART_MASK_FILES_PREFIX}.bed.gz
#' - ${PART_MASK_FILES_PREFIX}_masked.sync.gz

python3 ./snape-files/MaskSYNC_snape.py \
    --sync ${SYNC_FILE} \
    --output ${PART_MASK_FILES_PREFIX} \
    --indel ${MP_INFO_PREFIX}.indel \
    --coverage ${MP_INFO_PREFIX}.cov \
    --te ${REPEATS} \
    --mincov ${MIN_COV} \
    --maxcov ${MAX_COV} \
    --maxsnape 0 \
    --SNAPE

check_exit_status "MaskSYNC_snape (partial)" $?

# This removes redundant _masked ending to this file:
mv ${PART_MASK_FILES_PREFIX}_masked.sync.gz ${PART_MASK_FILES_PREFIX}.sync.gz




rm ${MP_INFO_PREFIX}* ${REPEATS}

rm -r snape-files


#' ========================================================================
#' Handle output files

mv ${MASK_FILES_PREFIX}.sync.gz ${TARGET}/

mkdir ${PART_MASK_FILES_PREFIX}
mv ${PART_MASK_FILES_PREFIX}.sync.gz ${MASK_FILES_PREFIX}.snp.gz \
    ./${PART_MASK_FILES_PREFIX}/
tar -czf ${PART_MASK_FILES_PREFIX}.tar.gz ${PART_MASK_FILES_PREFIX}
mv ${PART_MASK_FILES_PREFIX}.tar.gz ${TARGET}/
rm -r ${PART_MASK_FILES_PREFIX}

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

