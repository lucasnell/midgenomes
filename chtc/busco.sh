#!/bin/bash


# have job exit if any command returns with non-zero exit status (aka failure)
set -e

# This takes as input the file name from /staging/lnell/ you want to analyze
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi

cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

echo $GENOME
echo -e "\n\n"
OUTDIR=busco__${GENOME/.fasta/}

busco \
    -m genome \
    -l diptera_odb10 \
    -i ${GENOME} \
    -o ${OUTDIR} \
    --cpu 24

# ~~~~~~~~~~~~~
# For now, we'll just delete the main directory. Change this when you
# finalize the pipeline.
# ~~~~~~~~~~~~~
# tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
# mv ${OUTDIR}.tar.gz /staging/lnell/busco/
rm -r ${OUTDIR} ${GENOME}










# #### If doing this on your laptop, use this code:
#  (I don't recommend doing it on your laptop, though. It takes a while.)
#
# docker pull ezlabgva/busco:v5.2.2_cv2
#
#
# cd ~/_data/busco_wd
#
# docker run \
#     -v $(pwd):/busco_wd \
#     -u $(id -u) \
#     ezlabgva/busco:v5.2.2_cv2
#
#
# # Find the container id using
# docker ps
#
# # then open up the container with...
# docker exec -it CONTAINER_ID /bin/sh
#
# busco -m genome -l diptera_odb10 \
#     -i besst-scaffolds.fasta \
#     -o busco_out --cpu 6
#
#



