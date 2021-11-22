#!/bin/bash


# #### If doing this on your laptop, use this code:
#
# docker pull ezlabgva/busco:v5.2.2_cv1
#
# export wd=~/_data/busco_wd
#
# docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.2.2_cv1
#
#
# docker run \
#     -v $wd:/busco_wd \
#     ezlabgva/busco:v5.2.2_cv1 \
#     busco \
#         -m genome \
#         -l diptera_odb10 \
#         -i busco_wd/redundans_contigs.fasta \
#         -o busco_wd/busco_out \
#         --cpu 6
#
#



# For the cluster, use this:

for GENOME in contigs_shasta polished_hap1 polished_hap2
do
    cp /staging/lnell/${GENOME}.fasta.gz ./
    gunzip ${GENOME}.fasta.gz
    OUTDIR=busco__${GENOME}

    echo $GENOME
    echo -e "\n\n"

    busco \
        -m genome \
        -l diptera_odb10 \
        -i ${GENOME}.fasta \
        -o ${OUTDIR} \
        --cpu 24

    echo -e "\n\n\n\n"
    echo -e "============================================================"
    echo -e "============================================================"
    echo -e "\n\n\n\n"

    tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
    mv ${OUTDIR}.tar.gz /staging/lnell/
    rm -r ${OUTDIR} ${GENOME}.fasta
done



