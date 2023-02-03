#!/bin/bash

#' Convert a bunch of files from full species names to abbreviated names
#' with a description of what the file is.

# Description:
desc="contigs.fasta.gz"

fastas=$(find . -maxdepth 1 -type f ! -name "*_${desc}")

for f in $fastas; do
    g=$(python -c "x='$f'.split('_'); print(x[0][2] + x[1][0:5])")
    # Only switch commenting on below two lines once you've checked this
    echo mv $f ${g}_${desc}
    # mv $f ${g}_${desc}
done
