#!/bin/bash

export THREADS=32

# installation steps for Mambaforge
cp /staging/lnell/Mambaforge-Linux-x86_64.sh ./
export HOME=$PWD
export PATH
sh Mambaforge-Linux-x86_64.sh -b -p $PWD/mamba3
export PATH=$PWD/mamba3/bin:$PATH
rm Mambaforge-Linux-x86_64.sh

export CONDA_PREFIX=$PWD/mamba3

# Install snakemake:
mamba install -c conda-forge -c bioconda --yes snakemake
# To narrow FASTA files below
mamba install -c bioconda fastx_toolkit

# # example ------------------------------
# # wget https://github.com/a-ludi/dentist/releases/download/v3.0.0/dentist-example.tar.gz
# cp /staging/lnell/dentist-example.tar.gz ./
# tar -xzf dentist-example.tar.gz
# cd dentist-example
#
# snakemake --configfile=snakemake.yml --use-conda --cores=${THREADS}
# md5sum -c checksum.md5
# # --------------------------------------


# Argument from submit file dictates input FASTA and output name:

export OUTDIR=filled_dentist__$1
export SCAFFOLDS=$1.fasta
if [ ! -f /staging/lnell/${SCAFFOLDS}.gz ]; then
    echo "/staging/lnell/${SCAFFOLDS}.gz does not exist." 1>&2
    exit 1
fi

mkdir ${OUTDIR}
chmod +w -R ${OUTDIR}

# Move all files from submit node to OUTDIR:
mv dentist.v3.0.0.x86_64.tar.gz dentist.yml snakemake.yml ./${OUTDIR}/
cd ${OUTDIR}
# Adjust the snakemake.yml file for the input FASTA:
sed -i "s/INPUT_SCAFFOLDS_FILE/${SCAFFOLDS}/" snakemake.yml

# Now copy over scaffolds file:
cp /staging/lnell/${SCAFFOLDS}.gz ./ && gunzip ${SCAFFOLDS}.gz

# It's important that this is a FASTA file!
READS=basecalls_guppy-5.0.11.fasta
cp /staging/lnell/${READS}.gz ./ && gunzip ${READS}.gz
# Apparently the fasta2DAM > DPsplit step can't handle a file with periods
mv ${READS} basecalls_guppy.fasta
export READS=basecalls_guppy.fasta



export TMPDIR=$(pwd)/tmp
mkdir ${TMPDIR}


# Prepare dentist:
tar -xzf dentist.v3.0.0.x86_64.tar.gz
rm dentist.v3.0.0.x86_64.tar.gz

# copy other necessary files:
cp -r -t . \
    dentist.v3.0.0.x86_64/Snakefile \
    dentist.v3.0.0.x86_64/envs \
    dentist.v3.0.0.x86_64/scripts


# Make sure the FASTA files aren't too wide
mv ${SCAFFOLDS} ${SCAFFOLDS/.fasta/_orig.fasta}
mv ${READS} ${READS/.fasta/_orig.fasta}
fasta_formatter -i ${SCAFFOLDS/.fasta/_orig.fasta} -w 80 \
    -o ${SCAFFOLDS}
fasta_formatter -i ${READS/.fasta/_orig.fasta} -w 80 \
    -o ${READS}
rm ${SCAFFOLDS/.fasta/_orig.fasta} ${READS/.fasta/_orig.fasta}


# This section doesn't seem needed:
# # Now clean up sequence names
# # ... by removing everything after the first comma in the seq names:
# sed -r 's/\,.+//' ${SCAFFOLDS/.fasta/_narrow.fasta} > ${SCAFFOLDS}
# # ... by removing everything after the first space in the seq names:
# sed -r 's/\ .+//' ${READS/.fasta/_narrow.fasta} > ${READS}
# rm ${SCAFFOLDS/.fasta/_narrow.fasta} ${READS/.fasta/_narrow.fasta}


snakemake --configfile=snakemake.yml --use-conda \
    --conda-prefix ${CONDA_PREFIX} --cores=${THREADS}


rm -rf dentist.v3.0.0.x86_64 ${READS} ${SCAFFOLDS}

cd ..

tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

rm -rf ${TMPDIR} mamba3 ${OUTDIR}



