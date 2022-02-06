#!/bin/bash

# Used to summarize output.
# Usage: ./summarize.sh [THREADS] [OUT_FASTA] [OUT_DIR] [OUT_CSV]

export THREADS=$1
export OUT_FASTA=$2
export OUT_DIR=$3
export OUT_CSV=$4

if [ ! -f ${OUT_FASTA} ]; then
    echo -e "\n\nERROR: summarize.sh can't find FASTA (${OUT_FASTA})." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi


eval "$(conda shell.bash hook)"

conda activate main-env

# This outputs basics about scaffold sizes.
# First find it in PATH or in current directory.
# Resulting *.out file will be used below to make a CSV file of basic stats
# I'll also print this to stdout and save the file in this directory.
if [ $(command -v summ-scaffs.py) ]
then
    summ-scaffs.py ${OUT_FASTA} > scaff_summary.out
elif [ -f ./summ-scaffs.py ]
then
    ./summ-scaffs.py ${OUT_FASTA} > scaff_summary.out
else
    echo -e "\n\nERROR: summ-scaffs.py not found." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi

cat scaff_summary.out


# This outputs BUSCO scores
# Note that I'm doing the same thing with the stdout as with summ-scaffs
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS} > busco.out
conda deactivate

cat busco.out


# -------------
# Save summ-scaffs.py and BUSCO output into a simple CSV file.
# -------------
# Header:
echo -n "file,size,scaffs,N50,min,max,total_N," > ${OUT_CSV}
echo "BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n" >> ${OUT_CSV}
# output file name:
echo -n "${OUT_DIR}", >> ${OUT_CSV}
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


exit 0
