#!/bin/bash


#'
#' Find the genes common among all species and extract into combined files.
#'


#' All species in the phylogeny, in alphabetical order:
export ALL_SPECIES=(Anopheles_stephensi \
    Belgica_antarctica \
    Chironomus_riparius \
    Chironomus_tentans \
    Chironomus_tepperi \
    Clunio_marinus \
    Culicoides_sonorensis \
    Parochlus_steinenii \
    Polypedilum_pembai \
    Polypedilum_vanderplanki \
    Propsilocerus_akamusi \
    Tanytarsus_gracilentus)

#' I'm dropping Simulium_vittatum because it has way too many duplicate genes.
#     Simulium_vittatum \

mkdir working
cd working

for SPECIES in ${ALL_SPECIES[@]}; do
    ODB_DIR=${SPECIES}_odb
    tar -xf /staging/lnell/phylo/species-files/${ODB_DIR}.tar.gz -C ./
    status=$?
    unset ODB_DIR
    if ((status != 0)); then
        echo "failed on $SPECIES" 1>&2
        break
    fi
    unset status
done


# #' To view numbers of single-copy genes by species:
# for SPECIES in ${ALL_SPECIES[@]}; do
#     echo $SPECIES
#     ls -1 ${SPECIES}_odb | wc -l
#     echo -ne "\n"
# done


#'
#' Use some python code to find the common genes among all species and write
#' to "common_genes.txt" file.
#'
python3 << EOF
import glob
spp = [x.replace("_odb", "") for x in glob.glob("*_odb")]
genes = [[]] * len(spp)
for i in range(len(spp)):
    genes[i] = glob.glob(spp[i] + "_odb/*.faa")
    for j in range(len(genes[i])):
        x = genes[i][j]
        genes[i][j] = x.replace(spp[i] + "_odb/", "").replace(".faa", "")
common_genes = set(genes[0])
for s in genes[1:]:
    common_genes.intersection_update(s)
common_genes = list(common_genes)
# Sanity check:
for g in common_genes:
    for s in spp:
        z = glob.glob(s + "_odb/" + g + ".faa")
        if len(z) != 1:
            print("Gene '%s' found %i times for species '%s'" % (g, len(z), s))
# Write to file:
with open("common_genes.txt", "w") as file:
    for g in common_genes:
        file.write(g + "\n")
EOF


# #' To view number of common genes:
# cat common_genes.txt | wc -l


mkdir common_genes
# (brackets below are so that the progress bar shows properly)
{
total_genes=$(cat common_genes.txt | wc -l)
x=1

while read -r gene; do
    for sp in ${ALL_SPECIES[@]}; do
        cat ${sp}_odb/${gene}.faa \
            | sed "s/>.*/>${sp}__${gene}/g" \
            >> common_genes/${gene}.faa
    done
    perc=$(python3 -c "z=${x} / ${total_genes} * 100; print(f'{z:.2f}')")
    echo -ne "${perc}%\r"
    x=$((x + 1))
done < common_genes.txt
echo -ne "\n"
}

mv common_genes.txt ./common_genes/
tar -czf common_genes.tar.gz common_genes

mv common_genes.tar.gz /staging/lnell/phylo/

rm -r common_genes *_odb
