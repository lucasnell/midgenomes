#!/usr/bin/env python3
"""
Convert gSYNC files (optionally gzipped) into a pileup file.

NOTE: does not work with indels!

usage:
./sync2pileup.py \
    -o <output file name> \
    <input sync file as *.sync or *.sync.gz>

"""
import sys
import os.path
import gzip
import argparse
from datetime import datetime




def make_error(err_msg):
    sys.stderr.write("ERROR: " + err_msg + "\n")
    sys.exit(1)




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description = "Combine SNAPE info from sync files")
    parser.add_argument("-o", "--out", required = True, 
                        help = "Output file name")
    parser.add_argument("sync_in", help="input sync file")
    
    args = parser.parse_args()
    out = args.out
    sync_in = args.sync_in
    
    if not sync_in.endswith((".sync", ".sync.gz")):
        make_error("Strange suffix to input file name (must be .sync or " +
                   ".sync.gz). Exiting.")
    if not os.path.exists(sync_in):
        make_error(sync_in + " does not exist. Exiting.\n")
    if not out.endswith((".pileup", ".pileup.gz", ".txt", ".txt.gz")):
        make_error("Strange suffix to output file name (must be .pileup, " +
                   ".pileup.gz, .txt, or .txt.gz). Exiting.")
    out_dir_file = os.path.split(out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        make_error("output dir (" + out_dir_file[0] + ") not found")
    if os.path.exists(out):
        make_error(out + " already exists. Overwriting not allowed.\n")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Converting sync file to " + out)
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading/writing (" + time_str + ")")
    
    if sync_in.endswith(".gz"):
        in_file = gzip.open(sync_in, "rt")
    else:
        in_file = open(sync_in, "r")
    
    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")
    
    # Bases in the order they're listed in the sync file's 4th column.
    bases = "ATCGN"
    # Vector to map characters to 0-based indices for where the allele 
    # count should be in the sync file's 4th column.
    base_map = [-1 for x in range(255)]
    for i, x in enumerate(bases):
        base_map[ord(x)] = i
        base_map[ord(x.lower())] = i
    
    
    # this gets changed each iteration to have "." used for reference:
    tmp_bases = [x for x in bases]
    
    for i, line in enumerate(in_file):
        l_split = line.rstrip().split("\t")
        if len(l_split) < 4:
            make_error("only " + str(len(l_split)) + " columns in file " + 
                       sync_in + " line " + str(i))
        if l_split[3] == ".:.:.:.:.:.":
            continue
        new_line = "\t".join(l_split[:3])
        ref = l_split[2]
        ref_idx = base_map[ord(ref)]
        tmp_bases[ref_idx] = "."
        counts = [int(x) for x in l_split[3].split(":")]
        new_line += "\t"
        new_line += str(sum(counts))
        new_line += "\t"
        for j in range(len(tmp_bases)):
            new_line += tmp_bases[j] * counts[j]
        new_line += "\t"
        new_line += "~" * sum(counts)
        new_line += "\n"
        b = out_file.write(new_line)
        j += 1
        tmp_bases[ref_idx] = bases[ref_idx]
    
    in_file.close()
    out_file.close()
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

