#!/bin/bash

#' This script calls npstat to estimate population-genetic parameters from
#' masked sync files based on Pool-seq data.
#'
#' Parameters it estimates:
#' - nucleotide diversity
#' - Watterson’s theta
#' - Tajima’s D
#'



# export READ_BASE=Vik-19_S41
# export N_ADULTS=40
# export WIN_SIZE=10000

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
source /staging/lnell/helpers.sh


export READ_BASE=$1


#' The second input indicates the number of adult midges in this sample.
#'
export N_ADULTS=$2


#' The third indicates the size of sliding window in kb.
#'
export WIN_SIZE=${3}000


# To make sure that sequences are in the right order
export GENOME=tany_scaffolds.fasta





#' ========================================================================
#' Inputs
#' ========================================================================

export IN_DIR=${READ_BASE}_snape_part_masked
export SNP=${READ_BASE}_snape_masked.snp
export SYNC=${READ_BASE}_snape_part_masked.sync.gz



#' ========================================================================
#' Outputs
#' ========================================================================

# Where to send everything when done:
export TARGET=/staging/lnell/dna/npstat
# Final files / directories
export OUT_DIR=${READ_BASE}_npstat_w${WIN_SIZE}
export OUT_FILE=${OUT_DIR}.stat
# Intermediate file:
export PILEUP=${READ_BASE}_snape_part_masked.pileup


#' ========================================================================
#' Prep for downstream steps.
#' ========================================================================

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/dna/snape/${IN_DIR}.tar.gz ./
check_exit_status "cp main tar file" $?
tar -xzf ${IN_DIR}.tar.gz
cd ${IN_DIR}
gunzip ${SNP}.gz
mv ${SNP} ${SYNC} ../
cd ..
rm -r ${IN_DIR}.tar.gz ${IN_DIR}

cp /staging/lnell/sync2pileup.py ./
chmod +x sync2pileup.py
./sync2pileup.py -o ${PILEUP} ${SYNC}
check_exit_status "sync2pileup" $?
rm ${SYNC} sync2pileup.py

cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp ${GENOME}" $?

# Scaffold names in the same order as the reference.
# Used to organize output files since SNAPE analyzes by scaffold.
# Make sure scaffold names don't have spaces in them!
export SCAFF_NAMES=($(grep "^>" ${GENOME} | sed 's/>//g' | sed 's/\s.*$//'))
rm ${GENOME}



#' ========================================================================
#' Split pileup and SNP files by scaffold
#' ========================================================================
mkdir tmp
cd tmp

# Pileup file first:
# This does the splitting:
awk '{if (last != $1) close(last); print >> $1; last = $1}' ../${PILEUP}
# Rename scaffold files and verify they exist in the reference
for scaff in *; do
    if [[ ! " ${SCAFF_NAMES[*]} " =~ " ${scaff} " ]]; then
        echo "Scaffold ${scaff} does not exist in reference! " 1>&2
        exit 1
    fi
    mv $scaff ../${scaff}.pileup
done
rm ../${PILEUP}


# Now same for SNP file
awk '{if (last != $1) close(last); print $2 >> $1; last = $1}' ../${SNP}
for scaff in *; do
    if [[ ! " ${SCAFF_NAMES[*]} " =~ " ${scaff} " ]]; then
        echo "Scaffold ${scaff} does not exist in reference! " 1>&2
        exit 1
    fi
    mv $scaff ../${scaff}.snp
done
rm ../${SNP}


# Now go back to main directory
cd ..
rm -r tmp

# # If you want to check which files are present
# for scaff in ${SCAFF_NAMES[@]}; do
#     if [ ! -f ${scaff}.pileup ]; then
#         echo "${scaff}.pileup not found"
#     fi
#     if [ ! -f ${scaff}.snp ]; then
#         echo "${scaff}.snp not found"
#     fi
# done


#' ========================================================================
#' npstat
#' ========================================================================

# Shell script to run npstat on a single scaffold.
# Each one takes a while to run, so I want to do this in parallel.
cat << EOF > one_npstat.sh
#!/bin/bash
export scaff=\${1}
if [ ! -f \${scaff}.pileup ] || [ ! -f \${scaff}.snp ]; then
    rm \${scaff}.pileup \${scaff}.snp 2> /dev/null
    exit 1
fi
npstat -n $((N_ADULTS * 2)) -l ${WIN_SIZE} -maxcov 10000 -nolowfreq 0 \\
    -snpfile \${scaff}.snp \${scaff}.pileup \\
    &> \${scaff}_npstat.out
status=\$?
if [ "\$status" != "0" ]; then
    exit 2
fi
rm \${scaff}.pileup \${scaff}.snp \${scaff}_npstat.out

# Keep only columns that are useful given my input parameters
cut -f 1,2,4-10 \${scaff}.pileup.stats \\
    > \${scaff}.stats
rm \${scaff}.pileup.stats

exit 0
EOF

chmod +x one_npstat.sh


cat << EOF > mp_npstat.py
#!/usr/bin/env python3
import glob
import subprocess as sp
import multiprocessing as mp
import sys
import argparse


def work(pile_file):
    """Defines the work unit on an input file"""
    ret = sp.call(["./one_npstat.sh", pile_file])
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run npstat in parallel")
    parser.add_argument("inputs", nargs = "+", help="input scaffold names")
    args = parser.parse_args()
    tasks = args.inputs
    with mp.Pool(processes=${THREADS}) as pool:
        exit_codes = pool.map(work, tasks)
    not_found = [tasks[j] for j,i in enumerate(exit_codes) if i == 1]
    failed = [tasks[j] for j,i in enumerate(exit_codes) if i == 2]
    if len(not_found) > 0:
        sys.stderr.write("WARNING: these scaffolds were not found : {}\n".format(not_found))
    if len(failed) > 0:
        sys.stderr.write("ERROR: npstat failed for these scaffolds : {}\n".format(failed))
        sys.exit(1)
    sys.exit(0)

EOF

chmod +x mp_npstat.py

./mp_npstat.py ${SCAFF_NAMES[@]}
status=$?
for scaff in ${SCAFF_NAMES[@]}; do
    if [ -f ${scaff}_npstat.out ]; then
        echo -e "\n\nnpstat output for failed attempt at scaffold '${scaff}':" 1>&2
        cat ${scaff}_npstat.out 1>&2
        echo -e "\n-------------------------------------------\n" 1>&2
    fi
done
check_exit_status "mp_npstat" $status
rm mp_npstat.py one_npstat.sh



# Combine output files.
# Add header from one output file (all are the same) to the final output.
stats_outs=($(ls *.stats))
echo -ne "sequence\t" >> ${OUT_FILE}
head -n 1 ${stats_outs[0]} >> ${OUT_FILE}
# I'm using `SCAFF_NAMES` so that the combined file is in the same order
# as the reference.
for scaff in ${SCAFF_NAMES[@]}; do
    if [ -f ${scaff}.stats ]; then
        # remove header, add the scaffold column, then add this to main output file
        tail -n+2 ${scaff}.stats \
            | awk -v SCAFF_NAME="${scaff}" '{print SCAFF_NAME"\t"$0}' - \
            >> ${OUT_FILE}
        rm ${scaff}.stats
    fi
done




#' Modified by LAN on 2022-04-28 bc output from `npstat` had columns 9-10 and
#' 11-12 switched compared to `npstat` manual
# 1. window number,
# 2. number of bases covered by sequences,
# 3. number of bases covered and with known outgroup allele,
# 4. average read depth,
# 5. number of segregating sites S,
# 6. Watterson estimator of θ,
# 7. Tajima’s Π estimator of heterozygosity,
# 8. Tajima’s D,
# 9. variance of the number of segregating sites,
# 10. variance of the Watterson estimator,
# 11. unnormalized Fay and Wu’s H,
# 12. normalized Fay and Wu’s H,
# 13. divergence per base (from outgroup),
# 14. nonsynonimous polymorphisms,
# 15. synonimous polymorphisms,
# 16. nonsynonimous divergence,
# 17. synonimous divergence,
# 18. α (fraction of substitutions fixed by positive selection).





#' ========================================================================
#' Handle output file
#' ========================================================================


gzip ${OUT_FILE}
mv ${OUT_FILE}.gz ${TARGET}/

cd ..
rm -r ${OUT_DIR}

