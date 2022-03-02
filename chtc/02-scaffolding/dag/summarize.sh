#!/bin/bash

# Used to summarize output.

export THREADS=24

export PREFIX=$1
export IN_FASTA=${PREFIX}.fasta
export OUT_DIR=${PREFIX}_summ
export OUT_CSV=${PREFIX}.csv


if [ ! -f /staging/lnell/${IN_FASTA}.gz ]; then
    echo -e "\n\nERROR: /staging/lnell/${IN_FASTA}.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi


# If we've already summarized this in the main CSV file,
# this job stops with exit code 0
export GREP_EXIT=$(gunzip -c /staging/lnell/scaffolds_all.csv.gz | \
                   grep -q "^${PREFIX},"; \
                   echo $?)

case ${GREP_EXIT} in

  0)
    echo -e "\n\nMESSAGE: Output already found in /staging/lnell/scaffolds_all.csv.gz."
    echo -e "Exiting...\n"
    exit 0
    ;;

  1)
    echo -e "\n\nMESSAGE: Output not found. Continuing with script...\n"
    ;;

  *)
    echo -e "\n\nERROR: Error occurred when searching /staging/lnell/scaffolds_all.csv.gz."
    echo -e "Exiting...\n"
    exit 1
    ;;
esac


mkdir ${OUT_DIR}
cd ${OUT_DIR}


cp /staging/lnell/${IN_FASTA}.gz ./ && gunzip ${IN_FASTA}.gz


eval "$(conda shell.bash hook)"

conda activate main-env

# This outputs basics about scaffold sizes.
# Resulting *.out file will be used below to make a CSV file of basic stats
# I'll also print this to stdout and save the file in this directory.
summ-scaffs.py ${IN_FASTA} | \
    tee scaff_summary.out


# This outputs BUSCO scores
# Note that I'm doing the same thing with the stdout as with summ-scaffs
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${IN_FASTA} \
    -o busco \
    --cpu ${THREADS} | \
    tee busco.out
conda deactivate

# -------------
# Save summ-scaffs.py and BUSCO output into a CSV file.
# Also add this info to a file inside `/staging/lnell` that stores it for a
# bunch of files.
# -------------

# Header:
echo -n "file,size,scaffs,N50,min,max,total_N," > ${OUT_CSV}
echo "BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n" >> ${OUT_CSV}

# output file identifier:
echo -n "${PREFIX}", >> ${OUT_CSV}
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

tail -n 1 ${OUT_CSV} | gzip >> /staging/lnell/scaffolds_all.csv.gz

cat ${OUT_CSV}

cd ..

# uncomment these lines if you want to save this directory:
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/


rm -r ${OUT_DIR}
