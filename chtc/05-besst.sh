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
# conda pack -n besst-envv --ignore-missing-files
# chmod 644 besst-env.tar.gz
# ls -sh besst-env.tar.gz

# CHANGE THESE FOR DIFFERING PIPELINES:
# The first adjusts inputs, the second the name of outputs:
export GENOME=haploid_purge_dups.fasta
export OUT_SUFFIX=purge_dups


cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

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
export OUTDIR=scaffolds_besst_${OUT_SUFFIX}
export OUTFASTA=${OUTDIR}.fasta

# conda setup
ENVNAME=besst-env
ENVDIR=$ENVNAME
export PATH
mkdir ${ENVDIR}
cp /staging/lnell/${ENVNAME}.tar.gz ./
tar -xzf ${ENVNAME}.tar.gz -C ${ENVDIR}
. ${ENVDIR}/bin/activate
rm ${ENVNAME}.tar.gz

cd BESST_RNA

python Main.py 1 \
    -c ${WD}/${GENOME} \
    -f ${WD}/${BAM} \
    -o ${WD}/${OUTDIR} \
    -z 10000


# -e 10 uses more reads than default of 3
# -k 50 includes very small contigs
# -z 10000 effectively turns repeat detection off (this was recommended)




cd ${WD}

rm -r BESST_RNA ${BAM} ${BAM}.bai ${GENOME} ${ENVNAME}


# Move just the scaffolds to the staging directory:
cp ./${OUTDIR}/pass1/Scaffolds-pass1.fa ./
mv Scaffolds-pass1.fa ${OUTFASTA}
gzip ${OUTFASTA}
mv ${OUTFASTA}.gz /staging/lnell/

# practice run shows that it makes the `${OUTDIR}Statistics.txt` file
# outside of the output directory,
# so we want to make sure it's in the output dir, too
mv ${OUTDIR}Statistics.txt ./${OUTDIR}

tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/
rm -r ${OUTDIR}
