#!/usr/bin/env python3
"""
Summarize coverage from mpileup using a sliding window

usage:
    ./window-mpileup.py \
    -r <reference assembly> \
    -s <step size (default 50)> \
    -w <window size (default 100)> \
    <input mpileup as *.txt or *.txt.gz>

"""


import sys
import os.path
import gzip
import argparse
import numpy as np
import math
import statistics as stats
from datetime import datetime




def summarize_ref(reference):
    
    if reference.endswith(".gz"):
        ref_file = gzip.open(reference,"rt")
    else:
        ref_file = open(reference, "r")
    
    sizes = []
    names = []
    i = -1
    
    for line in ref_file:
        if line.startswith(">"):
            sizes.append(0)
            seq_name = line[1:]
            seq_name = seq_name.strip(" \t\n\r")
            names.append(seq_name)
            i += 1
        else:
            if len(sizes) == 0:
                sys.stderr.write("FASTA doesn't start with header. Exiting.\n")
                sys.exit(1)
            # we don't want to include newline at the end in our counts:
            llen = len(line) - 1
            sizes[i] += llen
    
    ref_file.close()
    
    return sizes, names



def one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names):
    """
    Summarize one window. Updates the output file directly.
    """
    
    if len(cov_ij) == 0:
        out_line = names[seq_i] + "\t"
        out_line += str(starts[seq_i][win_j]) + "\t"
        out_line += str(ends[seq_i][win_j]) + "\t"
        out_line += "0\t" + "0\t" + "0.0\t" + "0.0\n"
        b = out_file.write(out_line)
        return
    
    size_ij = np.int64(ends[seq_i][win_j] - starts[seq_i][win_j] + 1)
    if len(cov_ij) > size_ij:
        sys.stderr.write("ERROR: For window " + str(win_j+1) + 
                         " in sequence '" + names[seq_i] + 
                         "', the number of coverage values exceeds " +
                         "the window size!\n")
        sys.exit(1)
    
    n_zeros = size_ij - np.int64(len(cov_ij))
    for i in range(n_zeros):
        cov_ij.append(np.float64(0))
    
    min_ij = np.int64(min(cov_ij))
    max_ij = np.int64(max(cov_ij))
    mean_ij = stats.mean(cov_ij)
    median_ij = stats.median(cov_ij)
    
    out_line = names[seq_i] + "\t"
    out_line += str(starts[seq_i][win_j]) + "\t"
    out_line += str(ends[seq_i][win_j]) + "\t"
    out_line += str(min_ij) + "\t"
    out_line += str(max_ij) + "\t"
    out_line += str(mean_ij) + "\t"
    out_line += str(median_ij) + "\n"
    b = out_file.write(out_line)
    
    return


# For testing:
# in_mpileup = "Lys-19_S14_bwa_mpileup.txt.gz"
# reference = "tany_scaffolds.fasta.gz"
# step = 50
# window = 100



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "Summarize mpileup")
    parser.add_argument("-r", "--reference", required = True,
                        help = "reference assembly file name")
    parser.add_argument("-s", "--step", type=int, default=50,
                        help = "sliding window step size (default = 50)")
    parser.add_argument("-w", "--window", type=int, default=100,
                        help = "max coverage before binning (default = 100)")
    parser.add_argument("in_mpileup", help="input mpileup file name")
    
    args = parser.parse_args()
    in_mpileup = args.in_mpileup
    reference = args.reference
    step = args.step
    window = args.window
    
    
    # ---------------
    # Check files:
    # ---------------
    if not os.path.exists(in_mpileup):
        sys.stderr.write(in_mpileup + " not found\n")
        sys.exit(1)
    mp_suffs = (".txt.gz", ".txt")
    if not in_mpileup.endswith(mp_suffs):
        sys.stderr.write("Only *.txt or *.txt.gz extension allowed. Exiting.\n")
        sys.exit(1)
    fasta_suffs = (".fasta.gz", ".fasta", ".fa.gz", ".fa")
    if not reference.endswith(fasta_suffs):
        sys.stderr.write("Strange suffix to input FASTA file. Exiting.\n")
        sys.exit(1)
    out_fn = in_mpileup.replace(".txt", "_window.txt")
    if os.path.exists(out_fn):
        sys.stderr.write(out_fn + " already exists. Overwriting not allowed.\n")
        sys.exit(1)
    
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Sliding windows of coverage for " + in_mpileup + "\n")
    
    # --------------
    # Read reference file:
    # --------------
    sizes, names = summarize_ref(reference)
    num_seqs = len(sizes)
    
    # --------------
    # Make output objects:
    # --------------
    # vector of windows per sequence:
    win_vec = np.array([math.ceil((ss - window + 1) / step) for ss in sizes])
    
    # starting and ending points (inclusive) for each window:
    starts = [np.arange(0, ss-window+1, step) for ss in sizes]
    ends = [np.append(np.arange(window-1, ss-step, step), ss-1) for ss in sizes]
    
    if in_mpileup.endswith(".gz"):
        mp_file = gzip.open(in_mpileup,"rt")
    else:
        mp_file = open(in_mpileup, "r")
    if out_fn.endswith(".gz"):
        out_file = gzip.open(out_fn,"wt")
    else:
        out_file = open(out_fn, "w")
    
    out_line = "scaff\t" + "start\t" + "end\t" + "min\t" + "max\t"
    out_line += "mean\t" + "median\n"
    # Assigning to `b` so it doesn't print # bytes printed
    b = out_file.write(out_line)
    
    # indices for sequence and window
    seq_i = np.int64(0)
    win_j = np.int64(0)
    # to calc. median, we need to retain all coverage values within a window:
    cov_ij = []
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("starting processing files (" + time_str + ")")
    
    for line in mp_file:
        
        scaff, pos, ref, cov, read_base, qual = line.split("\t")
        pos = np.int64(pos)
        # mpileup uses 1-based indices
        pos -= 1
        cov = np.float64(cov)
        
        if scaff != names[seq_i]:
            """
            If we need to move to a new sequence, then summarize this window
            and, if applicable, any remaining in the sequence.
            Then do this for all windows in the sequences between this one
            and the one you need.
            """
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            cov_ij = []
            win_j += 1
            while scaff != names[seq_i]:
                while win_j < len(ends[seq_i]):
                    one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, 
                                  names)
                    win_j += 1
                seq_i += 1
                win_j = np.int64(0)
                if seq_i >= len(names):
                    sys.stderr.write("ERROR: assembly and mpileup sequence names " +
                                     "do not match or are not in the same order!")
                    sys.exit(1)
        
        if pos > ends[seq_i][win_j]:
            """
            If we need to change window within the sequence, then similarly
            summarize windows until you get to the one you want.
            """
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            cov_ij = []
            win_j += 1
            while pos > ends[seq_i][win_j]:
                one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, 
                                  names)
                win_j += 1
                if win_j >= len(ends[seq_i]):
                    sys.stderr.write("ERROR: position " + str(pos) + 
                                     " exceeds highest possible index for " + 
                                     "sequence '" + names[seq_i] + "' (" + 
                                     str(sizes[seq_i]-1)  + ")\n")
                    sys.exit(1)
        
        # Once you're at the sequence and window you need, then append this
        # coverage to `cov_ij`
        cov_ij.append(cov)
    
    # Summarize last window
    one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
    # Summarize any remaining empty windows:
    cov_ij = []
    win_j += 1
    while seq_i < len(names):
        while win_j < len(ends[seq_i]):
            one_summarize(out_file, cov_ij, seq_i, win_j, starts, ends, names)
            win_j += 1
        seq_i += 1
        win_j = np.int64(0)
    
    mp_file.close()
    out_file.close()
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished processing files (" + time_str + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

