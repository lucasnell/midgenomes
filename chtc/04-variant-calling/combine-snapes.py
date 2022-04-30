#!/usr/bin/env python3
"""
Combine SNAPE info from multiple gSYNC files (optionally gzipped) into one.
It assumes they're all in the same order, as they should be because this
format should include a line per bp in the reference.
The SNAPE info is in the 5th column, and when this isn't present, it's
replaced with ".:.:.:.:.:.:.".

To remove all missing data:
gunzip -c POOLS.sync.gz | grep -v ".:.:.:.:.:.:." | gzip > POOLS_noblanks.sync.gz

usage:
./combine-snapes.py \
    -o <output file name> \
    <input sync files as *.sync or *.sync.gz>

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
    
    parser = argparse.ArgumentParser(
        description = "Combine SNAPE info from sync files")
    parser.add_argument("-o", "--out", required = True, 
                        help = "Output file name")
    parser.add_argument("inputs", nargs = "+", help="input sync files")
    
    args = parser.parse_args()
    out = args.out
    inputs = args.inputs
    
    if not out.endswith((".sync", ".sync.gz")):
        make_error("Strange suffix to output file name. Exiting.")
    out_dir_file = os.path.split(out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        make_error("output dir (" + out_dir_file[0] + ") not found")
    if os.path.exists(out):
        make_error(out + " already exists. Overwriting not allowed.\n")

    if len(inputs) < 2:
        make_error("At least 2 input files required")
    
    gzipped = all(x.endswith(".sync.gz") for x in inputs)
    if not (gzipped or all(x.endswith(".sync") for x in inputs)):
        make_error("All inputs must have the same extension, and it must " +
                   "be *.sync.gz or *.sync")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Combining sync files into " + out)
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading/writing (" + time_str + ")")
    
    if gzipped:
        files = [gzip.open(x, "rt") for x in inputs]
    else:
        files = [open(x, "r") for x in inputs]
    
    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")
    
    j = 0
    for line in it.zip_longest(*files, fillvalue=""):
        l_split = line[0].rstrip().split("\t")
        if len(l_split) < 4:
            make_error("only " + str(len(l_split)) + " columns in file " + 
                       inputs[0] + " line " + str(j))
        zline = l_split[:4]
        if len(l_split) < 5: 
            zline[3] = ".:.:.:.:.:.:."
        else:
            zline[3] = l_split[4]
        for i in range(1, len(line)):
            l_split = line[i].split("\t")
            if len(l_split) < 4:
                make_error("only " + str(len(l_split)) + " columns in file " + 
                           inputs[i] + " line " + str(j))
            if len(l_split) < 5: 
                zline.append(".:.:.:.:.:.:.")
            else:
                zline.append(l_split[4].rstrip())
        new_line = "\t".join(zline) + "\n"
        b = out_file.write(new_line)
        j += 1
    
    for x in files:
        x.close()
    out_file.close()
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

