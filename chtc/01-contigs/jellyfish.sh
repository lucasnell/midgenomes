

# # First combine all NextPolish-corrected reads into one FASTA file:
#
# conda activate main-env
#
# cd ~/_data/next_correct
# tar -xzf next_correct.tar.gz
# cd ./01.seed_cns.sh.work
#
# for f in $(ls -d seed_cns*); do
#     cd $f
#     sed -r 's/\ /_/g' cns.fasta \
#         | seqtk seq -l 80 -A - \
#         >> ../next_corrected.fasta
#     cd ..
#     echo $f "finished"
# done
#
# conda deactivate
#
# gzip next_corrected.fasta
# mv next_corrected.fasta.gz ~/_data/


# Now analyze using jellyfish:

cd ~/_data
docker run -it --rm=true --platform linux/amd64 \
    -v $(pwd):/data \
    lucasnell/tany_genomics:v0.5.8 /bin/bash

mamba create -y -c bioconda -c conda-forge -n jelly-env kmer-jellyfish=2.3.0
conda activate jelly-env

jellyfish count -C -t 4 -m 21 -o tany_kc.jf \
     -s 100M \
    /data/trimmed_MyKS-19-B_S18_L002_R1_001.fastq \
    /data/trimmed_MyKS-19-B_S18_L002_R2_001.fastq

#      -s 100M \
#     <(zcat /data/next_corrected.fasta.gz)

jellyfish histo -t 4 tany_kc.jf > tany_kc.histo

mv tany_kc* /data/jellyfish/ill/
# mv tany_kc* /data/jellyfish/ont/

# exit
