#!/usr/bin/env python3
"""
FIlter REads from SOrted SUMmary information (fire-sosum)

This script filters reads from a FASTQ file for a given coverage of the 
Tanytarsus gracilentus genome based on a summary information file 
that's been sorted to have the same order as the reads in the FASTQ file.
The sorted summary file provides average read quality, read length, and
read names, among other things.

usage:
    ./fire-sosum.py <summary file> <coverage> <FASTQ/FASTA file>

"""


import sys
import os.path
import numpy as np
import pandas as pd
import gzip


if __name__ == "__main__":
    
    if len(sys.argv) > 4:
        print("too many arguments")
        sys.exit(1)
    if not os.path.exists(sys.argv[1]):
        print(str(sys.argv[1]) + " not found")
        sys.exit(1)
    if not os.path.exists(sys.argv[3]):
        print(str(sys.argv[3]) + " not found")
        sys.exit(1)

    df = pd.read_csv(sys.argv[1], sep = "\t")
    scores = df["mean_qscore_template"].to_numpy()
    lengths = df["sequence_length_template"].to_numpy()
    
    coverage = float(sys.argv[2])
    # (estimate genome size is fixed at 100 Mb)
    needed_total = 100e6 * coverage
    if lengths.sum() < needed_total:
        print("Coverage too high")
        sys.exit(1)
    
    """
    Find average quality score that works.
    I'm working down from 10 bc that's plenty since average quality is a pretty
    poor metric anyway.
    """
    qual = 10.0
    while lengths[scores >= qual].sum() < needed_total:
        qual -= 1.0
        if qual <= 0:
            print("no quality filter needed")
            break
    
    # Now let's find max length that gets us the coverage we need:
    slen = 1000
    # first take big jumps up
    while lengths[lengths >= slen].sum() > needed_total:
        slen += 2500
    slen -= 2500
    # then take smaller jumps to narrow it down:
    while lengths[lengths >= slen].sum() > needed_total:
        slen += 500
    slen -= 500
    
    print("minimum average Q score: " + str(qual))
    print("minimum length: " + str(slen))
    
    ql_filter = np.where((lengths >= slen) & (scores >= qual))
    # not sure why, but `np.where` wraps output in tuple
    ql_filter = ql_filter[0]

    # Now filter FASTQ/FASTA file:
    fast_fn = str(sys.argv[3])
    fastq_suffs = (".fastq.gz", ".fastq", ".fq.gz", ".fq")

    # Change just the file name (keep directory the same) for output
    out_dir_file = os.path.split(fast_fn)
    out_name = os.path.join(out_dir_file[0], "filtered_" + out_dir_file[1])


    # keep track of the next read number we want to output:
    rn_next = 0
    # keep track of read number we're on in the file we're reading:
    rn_idx = 0
    # whether to write a line
    do_write = False
    if fast_fn.endswith(fastq_suffs):
        """
        Since the file is a FASTQ file, then it's all about line numbers
        because there are 4 per read.
        This should be faster than searching every line.
        We can do this because the summary file is already sorted.
        """
        if fast_fn.endswith(".gz"):
            with gzip.open(fast_fn,"rt") as in_fastq:
                with gzip.open(out_name,"wt") as out_fastq:
                    for i, line in enumerate(in_fastq):
                        new_read = i % 4 == 0
                        if new_read and rn_next >= ql_filter.size:
                            break
                        if new_read:
                            do_write = (rn_idx == ql_filter.item(rn_next))
                            if do_write:
                                rn_next+=1
                            rn_idx+=1
                        if do_write:
                            out_fastq.write(line)
        else:
            with open(fast_fn, "r") as in_fastq:
                with open(out_name, "w") as out_fastq:
                    for i, line in enumerate(in_fastq):
                        new_read = i % 4 == 0
                        if new_read and rn_next >= ql_filter.size:
                            break
                        if new_read:
                            do_write = (rn_idx == ql_filter.item(rn_next))
                            if do_write:
                                rn_next+=1
                            rn_idx+=1
                        if do_write:
                            out_fastq.write(line)
    else:
        print("Strange suffix to FASTQ file. Exiting.")
        sys.exit(1)

    print("Filtered reads written to " + out_name)
    
    
    sys.exit(0)




