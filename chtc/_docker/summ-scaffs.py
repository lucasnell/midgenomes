#!/usr/bin/env python3
"""
Summarize scaffolds

usage:
    ./summ-scaffs.py <input FASTA>

"""


import sys
import os.path
import gzip
import argparse


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "Summarize scaffolds")
    parser.add_argument("in_fasta", help="input FASTA file name")
    
    args = parser.parse_args()
    
    # ---------------
    # Check files:
    # ---------------
    if not os.path.exists(args.in_fasta):
        sys.stderr.write(args.in_fasta + " not found\n")
        sys.exit(1)
    fasta_suffs = (".fasta.gz", ".fasta", ".fa.gz", ".fa")
    if not args.in_fasta.endswith(fasta_suffs):
        sys.stderr.write("Strange suffix to input FASTA file. Exiting.\n")
        sys.exit(1)
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Summarize scaffolds for " + args.in_fasta + "\n")

    if args.in_fasta.endswith(".gz"):
        fasta_file = gzip.open(args.in_fasta,"rt")
    else:
        fasta_file = open(args.in_fasta, "r")
        
    total_size = 0
    total_N = 0
    sizes = []
    i = -1
    
    for line in fasta_file:
        if line.startswith(">"):
            sizes.append(0)
            i += 1
        else:
            if len(sizes) == 0:
                sys.stderr.write("FASTA doesn't start with header. Exiting.\n")
                sys.exit(1)
            # we don't want to include newline at the end in our counts:
            llen = len(line) - 1
            sizes[i] += llen
            total_size += llen
            total_N += line.count('n')
            total_N += line.count('N')
    
    fasta_file.close()
    
    sizes.sort(reverse = True)
    
    
    n50_threshold = round(float(total_size) / 2.0)
    cum_sum = 0;
    i = 0;
    while i < len(sizes):
        cum_sum += sizes[i]
        if cum_sum >= n50_threshold:
            break
        i+=1
    n50 = sizes[i]
    
    print("size = " + str(total_size))
    print(str(len(sizes)) + " scaffolds")
    print("N50 = " + str(n50))
    print("min = " + str(sizes[-1]))
    print("max = " + str(sizes[0]))
    print("total N = " + str(total_N))
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)




