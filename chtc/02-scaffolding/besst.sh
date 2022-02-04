#!/bin/bash

# Use HISAT2 to align RNAseq reads to assembly.
# Then use BESST_RNA to scaffold assembly.

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

# Extra threads are for HISAT2 and BUSCO only
export THREADS=24

. /app/.bashrc

# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi

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


mkdir ${OUT_DIR}

cd ${OUT_DIR}

cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz




# ============================================================================
# ============================================================================
# HISAT2
# ============================================================================
# ============================================================================

conda activate main-env

export BAM=rna_hisat2.bam
export RNA1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA2=trimmed_TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz


hisat2-build ${GENOME} tany_hisat_idx
# (Pipe to samtools is to convert it to an aligned BAM file)
hisat2 -x tany_hisat_idx -1 ${RNA1} -2 ${RNA2} -k 3 -p ${THREADS} \
        --pen-noncansplice 1000000 --no-unal | \
    samtools sort -O bam - > ${BAM}
samtools index ${BAM} ${BAM}.bai

# (`tany_hisat_idx*` are index files from `hisat2-build`)
rm ${RNA1} ${RNA2} tany_hisat_idx*

conda deactivate



# ============================================================================
# ============================================================================
# BESST_RNA
# ============================================================================
# ============================================================================


conda activate besst-env


# The BESST_RNA scripts are a part of this Docker container and are in
# `/app/BESST_RNA`, so I'll move them here:
cp -rp /app/BESST_RNA ./

cd BESST_RNA

python Main.py 1 \
    -c ../${GENOME} \
    -f ../${BAM} \
    -o ../ \
    -z 10000

conda deactivate


cd ..
rm -r BESST_RNA ${GENOME}


if [ ! -f ./pass1/Scaffolds-pass1.fa ]; then
    echo "BESST_RNA failed to produce output." 1>&2
    echo " Saving BAM files to staging and quitting." 1>&2
    mv ${BAM} ${BAM}.bai /staging/lnell/
    cd ..
    rm -r ${OUT_DIR} *Statistics.txt
    exit 1
fi


# -e 10 uses more reads than default of 3
# -k 50 includes very small contigs
# -z 10000 effectively turns repeat detection off (this was recommended)

# It uses lowercase n instead of N, which I'm changing for consistency
# with other scaffolders:
sed -i -e '/^[^>]/s/n/N/g' ./pass1/Scaffolds-pass1.fa


# Move just the scaffolds to the staging directory:
cp ./pass1/Scaffolds-pass1.fa ./
mv Scaffolds-pass1.fa ${OUT_FASTA}
# Keep the uncompressed version for summaries below
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/

# practice run shows that it makes the `${OUT_DIR}Statistics.txt` file
# outside of the output directory,
# so we want to make sure it's in the output dir, too
cd ..
if [ -f *Statistics.txt ]; then
    mv *Statistics.txt ./${OUT_DIR}/
fi
cd ${OUT_DIR}



# This outputs basics about scaffold sizes:
summ-scaffs.py ${OUT_FASTA}


# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS}
conda deactivate


cd ..

# Now save the whole directory
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}
