#!/usr/bin/env python3
"""
FIlter REads using GUppy SUmmary file (fire-gusu)

This script filters reads from a FASTQ file to achieve a given coverage
of a genome.
Filtering is based on a summary information file like the one output from guppy.
This script filters based on inputs, then randomly chooses reads to achieve
the given coverage.
If the filters are too strict, it will loosen them until they can produce the
desired coverage.
It does this in increments of 10% for both read quality and length.
Output can be made reproducible by setting the `seed` argument.


usage:
    ./fire-gusu.py -s <summary file> -c <coverage> -g <genome size (Mb)> \
        -o <output FASTQ> <input FASTQ>

"""


import sys
import os.path
import numpy as np
import pandas as pd
import gzip
import argparse
from datetime import datetime



def sort_summary_file(summ_fn, fastq_fn):
    """
    Sort a summary file in the same order as the reads in a FASTQ file.
    This usually saves time versus searching for the read names in the 
    FASTQ file itself.
    """
    
    summ_df = pd.read_csv(summ_fn, sep = "\t")
    # filtering for passing filters bc if they didn't they aren't in the FASTQ
    summ_df = summ_df.loc[summ_df.passes_filtering]
    
    col_names = summ_df.columns.tolist()
    
    # (these will be changed later)
    fq_read_names = summ_df["read_id"].tolist()
    
    # set read name as index, for use in sorting later:
    summ_df = summ_df.set_index("read_id")
    
    # index in `fq_read_names`
    j = 0
    if fastq_fn.endswith(".gz"):
        in_fastq = gzip.open(fastq_fn,"rt")
    else:
        in_fastq = open(fastq_fn, "r")
    
    for i, line in enumerate(in_fastq):
        if i % 4 == 0:
            if j >= len(fq_read_names):
                sys.stderr.write("summary file is missing read names\n")
                sys.exit(1)
             # removes all after first space
            just_name = line.partition(" ")[0]
            #  `[1:]` removed starting @
            fq_read_names[j] = just_name[1:]
            j+=1
    
    in_fastq.close()
    
    sort_summ_df = summ_df.reindex(fq_read_names)
    # check for whether all read names match:
    if sort_summ_df["batch_id"].isnull().values.any():
        sys.stderr.write(str(sort_summ_df["batch_id"].isnull().sum()) + 
                         " read names don't match\n")
        sys.exit(1)
    
    if not all(sort_summ_df.index == fq_read_names):
        sys.stderr.write("sorted read names should match those "+
                         "from FASTQ file!\n")
        sys.exit(1)
    
    # Turning index back to column:
    sort_summ_df.reset_index(inplace = True)
    # Reordering to original column order:
    sort_summ_df = sort_summ_df[col_names]
    
    return sort_summ_df



def qual_len_filter(min_q, min_l, seq_needed, scores, lengths):
    """
    Return a vector of indices for reads that pass quality and length filters.
    No randomly choosing reads yet.
    """
    
    if lengths.sum() < seq_needed:
        sys.stderr.write("No filter needed for this depth of coverage.\n")
        sys.exit(1)
    
    # Reduce filtering if it's too strict:
    incr_q = min_q * 0.1
    incr_l = min_l * 0.1
    while lengths[(scores >= min_q) & (lengths >= min_l)].sum() < seq_needed:
        min_q -= incr_q
        min_l -= incr_l
        if min_q <= 0:
            break
    
    print("minimum average Q score: " + str(min_q))
    print("minimum length: " + str(min_l))
    
    ql_filter = np.where((lengths >= min_l) & (scores >= min_q))
    # not sure why, but `np.where` wraps output in tuple
    ql_filter = ql_filter[0]
    
    return min_q, min_l, ql_filter



def rnd_filter(lengths, ql_filter, seq_needed, rng):
    """
    Randomly choose reads until we get the desired coverage.
    """
    
    # no longer considering lengths not in filter:
    lengths = lengths[ql_filter]
    # Max allowed total sequencing output:
    max_seq_needed = seq_needed + np.mean(lengths)
    # No reason to continue if we're already within our range:
    if lengths.sum() <= max_seq_needed:
        return ql_filter
    
    # Shuffle indices, then remove samples until you have only the
    # sequencing needed:
    rnd_inds = np.arange(len(ql_filter))
    rng.shuffle(rnd_inds)
    # bp of sequencing in final output:
    out_seq = lengths[rnd_inds].sum()
    # number of reads in final output:
    out_nreads = len(rnd_inds)
    while out_seq > max_seq_needed:
        out_nreads -= 1
        if out_nreads <= 0:
            out_nreads = 1
            break
        # not subtracting one from index here bc we already have above:
        out_seq -= lengths[rnd_inds[out_nreads]]
    
    print("exact sequencing output: " + str(out_seq))
    
    if out_nreads < len(rnd_inds):
        rnd_inds = rnd_inds[:out_nreads] # first `out_nreads` reads
    
    # sort so ordering is still the same as in FASTQ
    rnd_inds.sort()
    
    rnd_ql_filter = ql_filter[rnd_inds.tolist()]
    
    return rnd_ql_filter








def write_filtered_fastq(in_fastq_fn, out_fastq_fn, ql_filter):
    """
    Based on a list of indices for reads to use, produce a new FASTQ file
    with only the filtered reads.
    """
    
    # to verify that we're outputting correct number of reads:
    nreads_out = 0
    # keep track of the next read number we want to output:
    rn_next = 0
    # keep track of read number we're on in the file we're reading:
    rn_idx = 0
    # whether to write a line
    do_write = False
    if in_fastq_fn.endswith(".gz"):
        in_fastq = gzip.open(in_fastq_fn,"rt")
    else:
        in_fastq = open(in_fastq_fn, "r")
    if out_fastq_fn.endswith(".gz"):
        out_fastq = gzip.open(out_fastq_fn,"wt")
    else:
        out_fastq = open(out_fastq_fn, "w")
    
    for i, line in enumerate(in_fastq):
        new_read = i % 4 == 0
        if new_read and rn_next >= ql_filter.size:
            break
        if new_read:
            do_write = (rn_idx == ql_filter.item(rn_next))
            if do_write:
                nreads_out+=1
                rn_next+=1
            rn_idx+=1
        if do_write:
            out_fastq.write(line)
    
    in_fastq.close()
    out_fastq.close()
    
    print("Filtered reads written to " + out_fastq_fn)
    
    if nreads_out != len(ql_filter):
        print("WARNING: " + str(len(ql_filter)) + " reads requested but only "+
              str(nreads_out) + " output to new file.")
    
    return















if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "FIlter REads using GUppy " +
                                                   "SUmmary file (fire-gusu)")
    parser.add_argument("-s", "--summary", required = True,
                        help = "path to tab-delimited summary file")
    parser.add_argument("-c", "--coverage", required = True, type=float, 
                        help = "coverage of genome desired")
    parser.add_argument("-g", "--genome_size", required = True, type=float,
                        help = "size of genome in Mb")
    parser.add_argument("-q", "--quality", required = True, type=float, 
                        help = "minimum average read quality")
    parser.add_argument("-l", "--length", required = True, type=float,
                        help = "minimum read length")
    parser.add_argument("--seed", required = False, type=int,
                        help = "set seed for random number generator for " + 
                        "reproducible output")
    parser.add_argument("-o", "--out_fastq", required = False,
                        help = "output file name (defaults to "+
                        "filtered_<input name>)")
    parser.add_argument("in_fastq", help="input file name")
    
    args = parser.parse_args()
    
    # ---------------
    # Check files:
    # ---------------
    if not os.path.exists(args.summary):
        sys.stderr.write(args.summary + " not found\n")
        sys.exit(1)
    if not os.path.exists(args.in_fastq):
        sys.stderr.write(args.in_fastq + " not found\n")
        sys.exit(1)
    if os.path.exists(args.out_fastq):
        sys.stderr.write("This script refuses to overwrite files. Good day.\n")
        sys.exit(1)
    
    fastq_suffs = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
    if not args.in_fastq.endswith(fastq_suffs):
        sys.stderr.write("Strange suffix to input FASTQ file. Exiting.\n")
        sys.exit(1)
    if args.out_fastq:
        if not args.out_fastq.endswith(fastq_suffs):
            sys.stderr.write("Strange suffix to output FASTQ file. Exiting.\n")
            sys.exit(1)
        out_dir_file = os.path.split(args.out_fastq)
        if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
            sys.stderr.write("output dir (" + out_dir_file[0] + ") not found\n")
            sys.exit(1)
        out_fastq = args.out_fastq
    else:
        # Change just the file name (keep directory the same) for output
        out_dir_file = os.path.split(args.in_fastq)
        out_fastq = os.path.join(out_dir_file[0], "filtered_" + out_dir_file[1])
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("sorting summary file (" + time_str + ")")
    # summary file sorted in the same order as the input FASTQ
    sorted_summ_df = sort_summary_file(args.summary, args.in_fastq)
    
    scores = sorted_summ_df["mean_qscore_template"].to_numpy()
    lengths = sorted_summ_df["sequence_length_template"].to_numpy()
    # Total sequencing needed:
    seq_needed = args.genome_size * 1e6 * args.coverage
    
    # filter for mean quality and length.
    # function also uses the min quality and length actually used, since 
    # too-strict filters are softened here so that we always have the 
    # coverage we want
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("filtering by quality and length (" + time_str + ")")
    min_q, min_l, ql_filter = qual_len_filter(args.quality, args.length, 
                                              seq_needed, scores, lengths)
    
    # Randomly sample from ql_filter to get desired coverage:
    if args.seed:
        rng = np.random.default_rng(args.seed)
    else:
        rng = np.random.default_rng()
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("randomly filtering for coverage (" + time_str + ")")
    rnd_ql_filter = rnd_filter(lengths, ql_filter, seq_needed, rng)

    # Do the filtering and output new FASTQ file:
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("writing output file (" + time_str + ")")
    write_filtered_fastq(args.in_fastq, out_fastq, rnd_ql_filter)

    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("finished (" + time_str + ")")
    
    sys.exit(0)




