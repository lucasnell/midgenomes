#!/bin/bash

#'
#' Use BUSCO to find single-copy genes from OrthoDB (diptera_odb10) within each
#' assembly, then find the genes shared amongst all species.
#'


. /app/.bashrc

export THREADS=$(count_threads)

export INPUT_LOC=/staging/lnell/assemblies
export OUTPUT_LOC=/staging/lnell/phylo/odb
export SHARED_OUT_DIR=shared_odb


#' Abbreviations of all species in the phylogeny, in alphabetical order:
export ALL_SPECIES=(Aaegyp Asteph Bantar Cripar Ctenta Cmarin Cquinq Csonor
                    Mdomes Pstein Ppemba Pvande Pakamu Tgraci)


mkdir ${SHARED_OUT_DIR}
cd ${SHARED_OUT_DIR}




#' ====================================================================
#' ====================================================================
#' Run BUSCO to match assemblies to diptera_odb10
#' ====================================================================
#' ====================================================================


conda activate busco-env

busco --download diptera_odb10 > busco-download.out

export BUSCO_DOWNLOADS=$(pwd)/busco_downloads/lineages/diptera_odb10

echo -e "species\tcount" > single-copy-gene-counts.tsv
for ((i=0; i<${#ALL_SPECIES[@]}; i++)); do

    sp=${ALL_SPECIES[$i]}

    GENOME=${sp}_assembly.fasta
    ODB_DIR=${sp}_odb
    BUSCO_DIR=${sp}_busco

    cp ${INPUT_LOC}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
    check_exit_status "cp ${GENOME}" $?

    busco -m genome -i ${GENOME} -o ${BUSCO_DIR} \
        --cpu ${THREADS} -l ${BUSCO_DOWNLOADS} \
    1> ${sp}_busco.stdout \
    2> ${sp}_busco.stderr
    check_exit_status "busco ${sp}" $?

    cp -r ${BUSCO_DIR}/run_diptera_odb10/busco_sequences/single_copy_busco_sequences ./ \
        && mv single_copy_busco_sequences ${ODB_DIR} \
        && mv ${sp}_busco.stdout ${sp}_busco.stderr ./${ODB_DIR}/ \
        && tar -czf ${ODB_DIR}.tar.gz ${ODB_DIR} \
        && mv ${ODB_DIR}.tar.gz ${OUTPUT_LOC}/
    check_exit_status "handle output ${sp}" $?

    rm -r ${GENOME} ${BUSCO_DIR}

    #' Number of single-copy genes for each species:
    COUNT=$(ls -1q ${sp}_odb/*.faa | wc -l)
    echo -e "${sp}\t${COUNT}" >> single-copy-gene-counts.tsv

    unset -v GENOME ODB_DIR BUSCO_DIR COUNT sp

    #' Provide progess if in an interactive session.
    if [[ $- == *i* ]]; then
        echo "finished species $(( $i + 1 )) of ${#ALL_SPECIES[@]}"
    fi


done


conda deactivate



#' ====================================================================
#' ====================================================================
#' Find and extract shared genes
#' ====================================================================
#' ====================================================================

#'
#' Use some python code to find the shared genes among all species and write
#' to "shared-genes.txt" file.
#'
python3 << EOF
import glob

spp = "${ALL_SPECIES[@]}".split(" ")
genes = [[]] * len(spp)
for i in range(len(spp)):
    genes[i] = glob.glob(spp[i] + "_odb/*.faa")
    for j in range(len(genes[i])):
        x = genes[i][j]
        genes[i][j] = x.replace(spp[i] + "_odb/", "").replace(".faa", "")
shared_genes = set(genes[0])
for s in genes[1:]:
    shared_genes.intersection_update(s)
shared_genes = list(shared_genes)
# Sanity check:
for g in shared_genes:
    for s in spp:
        z = glob.glob(s + "_odb/" + g + ".faa")
        if len(z) != 1:
            print("Gene '%s' found %i times for species '%s'" % (g, len(z), s))
# Write to file:
with open("shared-genes.txt", "w") as file:
    for g in shared_genes:
        file.write(g + "\n")
EOF


total_genes=$(cat shared-genes.txt | wc -l)
echo "number of shared genes = ${total_genes}"



#'
#' Extract shared genes.
#' (Uses progress bar in an interactive session, and brackets below are for the
#'  bar to display properly.)
#'
{
x=1
while read -r gene; do
    for sp in ${ALL_SPECIES[@]}; do
        cat ${sp}_odb/${gene}.faa \
            | sed "s/>.*/>${sp}__${gene}/g" \
            >> ${gene}.faa
    done
    if [[ $- == *i* ]]; then
        perc=$(python3 -c "z=${x} / ${total_genes} * 100; print(f'{z:.2f}')")
        echo -ne "${perc}%\r"
        x=$((x + 1))
    fi
done < shared-genes.txt
if [[ $- == *i* ]]; then echo -ne "\n"; fi
}
#' These are no longer needed:
rm -r shared-genes.txt ${ALL_SPECIES[@]/%/_odb} busco_downloads



#' ====================================================================
#' ====================================================================
#' Handle output
#' ====================================================================
#' ====================================================================


cd ..

tar -czf ${SHARED_OUT_DIR}.tar.gz ${SHARED_OUT_DIR} \
    && mv ${SHARED_OUT_DIR}.tar.gz ${OUTPUT_LOC}/
check_exit_status "handle shared genes output" $?

rm -r ${SHARED_OUT_DIR}

