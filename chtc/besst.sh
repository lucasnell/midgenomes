#!/bin/bash

# # Creating the conda environment is a bit different here bc I have to use
# # old versions for everything:
# conda create --name besst-env -c free python=2.7.10
# conda activate besst-env
# conda config --add channels r
# conda config --add channels bioconda
# conda config --add channels free
# conda install pysam matplotlib
# # This one requires previous version:
# conda install networkx==1.10
# conda deactivate
# conda install -c conda-forge conda-pack
# conda pack -n besst-env --ignore-missing-files
# chmod 644 besst-env.tar.gz
# ls -sh besst-env.tar.gz

# Extra threads are for BUSCO only
export THREADS=4

# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


# The input file dictates the output name:
export OUT_SUFFIX=besst

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME/.fasta/}_${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta



# The genome file affects which BAM file to use:
export BAM=rna_hisat2__${GENOME/.fasta/}.bam
if [ ! -f /staging/lnell/${BAM} ]; then
    echo "/staging/lnell/${BAM} does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${BAM} ./
cp /staging/lnell/${BAM}.bai ./


tar -xzf BESST_RNA.tar.gz
rm BESST_RNA.tar.gz

export WD=$(pwd)

# # conda setup
# ENVNAME=besst-env
# ENVDIR=$ENVNAME
# export PATH
# mkdir ${ENVDIR}
# cp /staging/lnell/${ENVNAME}.tar.gz ./
# tar -xzf ${ENVNAME}.tar.gz -C ${ENVDIR}
# . ${ENVDIR}/bin/activate
# rm ${ENVNAME}.tar.gz

conda activate besst-env

cd BESST_RNA

python Main.py 1 \
    -c ${WD}/${GENOME} \
    -f ${WD}/${BAM} \
    -o ${WD}/${OUT_DIR} \
    -z 10000

conda deactivate


cd ..
rm -r BESST_RNA ${BAM} ${BAM}.bai ${GENOME}


if [ ! -f ./${OUT_DIR}/pass1/Scaffolds-pass1.fa ]; then
    echo "BESST_RNA failed to produce output." 1>&2
    rm -r ${OUT_DIR} ${OUT_DIR}Statistics.txt
    exit 1
fi


# -e 10 uses more reads than default of 3
# -k 50 includes very small contigs
# -z 10000 effectively turns repeat detection off (this was recommended)

# It uses lowercase n instead of N, which I'm changing for consistency
# with other scaffolders:
cd ./${OUT_DIR}/pass1
sed -i -e '/^[^>]/s/n/N/g' Scaffolds-pass1.fa


cd ..
cd ..


# Move just the scaffolds to the staging directory:
cp ./${OUT_DIR}/pass1/Scaffolds-pass1.fa ./
mv Scaffolds-pass1.fa ${OUT_FASTA}
# Keep the uncompressed version for summaries below
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/

# practice run shows that it makes the `${OUT_DIR}Statistics.txt` file
# outside of the output directory,
# so we want to make sure it's in the output dir, too
mv ${OUT_DIR}Statistics.txt ./${OUT_DIR}


# This outputs basics about scaffold sizes:
summ_scaffs ${OUT_FASTA}

export BUSCO_OUT=busco_${OUT_DIR}


# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o ${BUSCO_OUT} \
    --cpu ${THREADS}
conda deactivate


# ~~~~~~~~~~~~~
# For now, we'll just delete the main and BUSCO directories.
# Change this when you finalize the pipeline.
# ~~~~~~~~~~~~~
# # Now the whole directories
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/
# tar -czf ${BUSCO_OUT}.tar.gz ${BUSCO_OUT}
# mv ${BUSCO_OUT}.tar.gz /staging/lnell/
rm -r ${OUT_DIR} ${BUSCO_OUT} ${OUT_FASTA}
