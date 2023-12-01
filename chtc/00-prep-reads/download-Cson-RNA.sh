#!/bin/bash

# Download RNAseq data for Culicoides sonorensis from SRA
# This does not need to be run on the `midgenomes` docker container.
# It should be run in a vanilla CentOS environment.
# (For midgenomes docker, use
#  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz)

mkdir working
cd working


export accessions=(ERR637904 ERR637905 ERR637906 ERR637907 ERR637908 ERR637909 \
                   ERR637910 ERR637911 ERR637912 ERR637913 ERR637914)

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64.tar.gz
tar -xzf sratoolkit.3.0.0-centos_linux64.tar.gz
rm sratoolkit.3.0.0-centos_linux64.tar.gz
export PATH=${PATH}:$(pwd)/sratoolkit.3.0.0-centos_linux64/bin

vdb-config --interactive
# just press x

# break loop if one fails:
check_exit_status () {
    if [ "$2" != "0" ]; then
        echo "Step $1 failed with exit status $2" 1>&2
        break 2> /dev/null
        return 0
    fi
    echo "Checked step $1"
    return 0
}



for acc in ${accessions[@]}; do
    prefetch --type fastq --progress $acc
    check_exit_status "downloading $acc" $?
done

export wd=$(pwd)

for acc in ${accessions[@]}; do
    cd ${acc}
    tar -cf ${acc}.tar *.fastq.gz
    check_exit_status "combining $acc" $?
    cd ${wd}
done


export destination=/staging/lnell/Cson-RNA
mkdir ${destination}/

for acc in ${accessions[@]}; do
    cd ${acc}
    mv ${acc}.tar ${destination}/
    check_exit_status "moving $acc" $?
    cd ${wd}
done


cd ..
rm -r working

