#!/usr/bin/env python3
"""
Combine multiple gSYNC files (optionally gzipped) into one.
It assumes they're all in the same order, as they should be because this
format should include a line per bp in the reference.
It also removes any instances of a 5th column in the output, which 
can happen when converting from SNAPE-pooled output to gSYNC.

For use in `poolfstat`, you should also run this to remove any missing data:
gunzip -c POOLS.sync.gz | grep -v ".:.:.:.:.:." | gzip > POOLS_noblanks.sync.gz

usage:
./combine-syncs.py -o <output file name> <input sync files as *.sync or *.sync.gz>

"""
import sys
import os.path
import gzip
import argparse
import itertools as it
from datetime import datetime



def make_error(err_msg):
    sys.stderr.write("ERROR: " + err_msg + "\n")
    sys.exit(1)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "Combine sync files")
    parser.add_argument("-o", "--out", required = True, 
                        help = "output file name")
    parser.add_argument("inputs", nargs = "+", help="input sync files")
    
    args = parser.parse_args()
    
    if not args.out.endswith((".sync", ".sync.gz")):
        make_error("Strange suffix to output file name. Exiting.")
    out_dir_file = os.path.split(args.out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        make_error("output dir (" + out_dir_file[0] + ") not found")
    if os.path.exists(args.out):
        make_error(args.out + " already exists. Overwriting not allowed.\n")
    out = args.out
    
    inputs = args.inputs
    if len(inputs) < 2:
        make_error("At least 2 input files required")
    
    gzipped = all(x.endswith(".sync.gz") for x in inputs)
    
    if not (gzipped or all(x.endswith(".sync") for x in inputs)):
        make_error("all inputs must have *.sync.gz or *.sync extensions")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Combining sync files into " + out + "\n")
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("opening files (" + time_str + ")")
    
    if gzipped:
        files = [gzip.open(x, "rt") for x in inputs]
    else:
        files = [open(x, "r") for x in inputs]
    
    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading/writing (" + time_str + ")")
    
    for line in it.zip_longest(*files, fillvalue=""):
        zline = line[0].rstrip().split("\t")[:4]
        for x in line[1:]:
            zline.append(x.split("\t")[3].rstrip())
        new_line = "\t".join(zline) + "\n"
        b = out_file.write(new_line)
    
    for x in files:
        x.close()
    out_file.close()
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished reading/writing (" + time_str + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

