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

gunzip redundans_contigs.fasta.gz

busco \
    -m genome \
    -l diptera_odb10 \
    -i redundans_contigs.fasta \
    -o busco_out \
    --cpu 12


tar -czf busco_out.tar.gz busco_out
rm -r busco_out redundans_contigs.fasta
