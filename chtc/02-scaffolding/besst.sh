#!/bin/bash

# Use HISAT2 to align RNAseq reads to assembly.
# Then use BESST_RNA to scaffold assembly.

# Extra threads are for HISAT2 only
export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
source /staging/lnell/helpers.sh


# Where to send and receive files from/to
export TARGET=/staging/lnell/assemblies

# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f ${TARGET}/${GENOME}.gz ]; then
    echo -e "\n\nERROR: ${TARGET}/${GENOME}.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
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

# # If the output FASTA already exists, this job stops with exit code 0
# if [ -f ${TARGET}/${OUT_FASTA}.gz ]; then
#     echo -e "\n\nMESSAGE: ${TARGET}/${OUT_FASTA}.gz already exists."
#     echo -e "Exiting...\n"
#     exit 0
# fi


mkdir ${OUT_DIR}

cd ${OUT_DIR}

cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz
check_exit_status "cp genome" $?



# ============================================================================
# ============================================================================
# HISAT2
# ============================================================================
# ============================================================================

conda activate main-env

export BAM=rna_hisat2.bam
export RNA1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA2=trimmed_TanyAdult_S1_L002_R2_001.fastq
if [ ! -f /staging/lnell/rna/trimmed_TanyAdult_S1.tar ]; then
    echo -e "\n\nERROR: /staging/lnell/rna/trimmed_TanyAdult_S1.tar does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi

cp /staging/lnell/rna/trimmed_TanyAdult_S1.tar ./ \
    && tar -xf trimmed_TanyAdult_S1.tar \
    && rm trimmed_TanyAdult_S1.tar
check_exit_status "move rna tar" $?

gunzip ${RNA1}.gz && gunzip ${RNA2}.gz
check_exit_status "gunzip rna" $?

hisat2-build ${GENOME} tany_hisat_idx
check_exit_status "hisat2-build" $?

# (Pipe to samtools is to convert it to an aligned BAM file)
hisat2 -x tany_hisat_idx -1 ${RNA1} -2 ${RNA2} -k 3 -p ${THREADS} \
        --pen-noncansplice 1000000 --no-unal | \
    samtools sort -O bam - > ${BAM}
check_exit_status "hisat2" $?
samtools index ${BAM} ${BAM}.bai
check_exit_status "samtools index" $?


# (`tany_hisat_idx*` are index files from `hisat2-build`)
rm ${RNA1} ${RNA2} tany_hisat_idx*



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

# -z 10000 effectively turns repeat detection off (this was recommended)
# -k 100 includes contigs >= 100 bp

python Main.py 1 \
    -c ../${GENOME} \
    -f ../${BAM} \
    -o ../ \
    -z 10000 \
    -k 100
check_exit_status "BESST_RNA" $?

conda deactivate


cd ..
rm -r BESST_RNA ${GENOME}


if [ ! -f ./pass1/Scaffolds-pass1.fa ]; then
    echo -e "\n\nERROR: BESST_RNA failed to produce output." 1>&2
    echo -e "Exiting...\n" 1>&2
    cd ..
    rm -r ${OUT_DIR} *Statistics.txt
    exit 1
fi




# It uses lowercase n instead of N, which I'm changing for consistency
# with other scaffolders:
sed -i -e '/^[^>]/s/n/N/g' ./pass1/Scaffolds-pass1.fa


# Move just the scaffolds to the staging directory:
cp ./pass1/Scaffolds-pass1.fa ./
mv Scaffolds-pass1.fa ${OUT_FASTA}
# Keep the uncompressed version for output directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

# practice run shows that it can make the `${OUT_DIR}Statistics.txt` file
# outside of the output directory,
# so we want to make sure it's in the output dir, too
cd ..
if [ -f *Statistics.txt ]; then
    mv *Statistics.txt ./${OUT_DIR}/
fi
cd ${OUT_DIR}

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_FASTA/.fasta/} | \
    tee ${OUT_FASTA/.fasta/}.csv


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}
