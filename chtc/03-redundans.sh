
export indir=~/_data/nanopore_assembly
export outdir=~/_data

docker run \
    -v $indir:/inputs:rw \
    -v $outdir:/output:rw \
    -it lpryszcz/redundans \
    /root/src/redundans/redundans.py \
        -v -t 4 \
        -l inputs/arives.guppy2.3.5.fastq.gz \
        -f inputs/midge.contigs.fasta \
        -o output/redundans

