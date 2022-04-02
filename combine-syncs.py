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
./combine-syncs.py \
    -o <output file name> \
    -t <# threads> \
    [--lowmem] \
    <input sync files as *.sync or *.sync.gz>

"""
import sys
import os.path
import gzip
import argparse
import itertools as it
from datetime import datetime
import multiprocessing as mp




def make_error(err_msg):
    sys.stderr.write("ERROR: " + err_msg + "\n")
    sys.exit(1)



def read_file_to_vec(i, file_name):
    """
    Read a single file's info into a vector.
    Arg `i` indicates the index for which file this is in the vector of inputs.
    It's only used to check whether `i == 0`, in which case the first 4 
    columns are kept.
    If `i != 0`, then just the 4th column is kept.
    """
    if file_name.endswith(".sync.gz"):
        file = gzip.open(file_name, "rt")
    elif file_name.endswith(".sync"):
        file = open(file_name, "r")
    else:
        make_error("file " + file_name +
                   " does not have *.sync.gz or *.sync extension")
    
    out_vec = []
    
    if i == 0:
        for j, line in enumerate(file):
            split_line = line.rstrip().split("\t")
            if len(split_line) < 4:
                make_error("only " + str(len(split_line)) + 
                           " columns in file " + file_name + " line " + str(j))
            split_line = split_line[:4]
            out_vec.append("\t".join(split_line))
    else:
        for j, line in enumerate(file):
            split_line = line.split("\t")
            if len(split_line) < 4:
                make_error("only " + str(len(split_line)) + 
                           " columns in file " + file_name + " line " + str(j))
            out_vec.append(split_line[3].rstrip())
    
    file.close()
    
    return out_vec





def lowmem_read_write(inputs, out):
    
    """
    Read and write simultaneously using the low-memory option for this script.
    """
    
    gzipped = all(x.endswith(".sync.gz") for x in inputs)
    if not (gzipped or all(x.endswith(".sync") for x in inputs)):
        make_error("For the low-memory option, all inputs must have the " +
                   "same extension, and it must be *.sync.gz or *.sync")
    
    print("  (using low-memory option)\n")
    
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
        for i in range(1, len(line)):
            l_split = line[i].split("\t")
            if len(l_split) < 4:
                make_error("only " + str(len(l_split)) + " columns in file " + 
                           inputs[0] + " line " + str(j))
            zline.append(l_split[3].rstrip())
        new_line = "\t".join(zline) + "\n"
        b = out_file.write(new_line)
        j += 1
    
    for x in files:
        x.close()
    out_file.close()







def highmem_read_write(inputs, out, threads):
    
    """
    Read inputs files into memory (potentially using >1 threads), 
    then write them to the output file.
    """
    
    print("")
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started reading (" + time_str + ")")
    
    pool = mp.Pool(threads)
    results = [pool.apply(read_file_to_vec, args=(i, x)) for i, x in enumerate(inputs)]
    pool.close()
    
    # Check that results all have same length:
    n_rows = len(results[0])
    if not all([len(x) == n_rows for x in results[1:]]):
        lens = [len(x) for x in results]
        lens_names = [x + " = " + str(y) for x, y in zip(inputs, lens)]
        make_error("Input files must have the same number of rows. Here " +
                   "are the row numbers for your input files:\n" +
                   "\n".join(lens_names))
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("started writing (" + time_str + ")")
    
    if out.endswith(".gz"):
        out_file = gzip.open(out, "wt")
    else:
        out_file = open(out, "w")
    
    n_files = len(results)
    for i in range(n_rows):
        lines_i = (results[j][i] for j in range(n_files))
        new_line = "\t".join(lines_i) + "\n"
        b = out_file.write(new_line)
    
    out_file.close()











if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "Combine sync files")
    parser.add_argument("-o", "--out", required = True, 
                        help = "Output file name")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help = "Number of threads")
    parser.add_argument("--lowmem", action="store_true",
                        help = str("Use the low-memory option. " +
                        "Not compatible with >1 threads."))
    parser.add_argument("inputs", nargs = "+", help="input sync files")
    
    args = parser.parse_args()
    out = args.out
    threads = args.threads
    lowmem = args.lowmem
    inputs = args.inputs
    
    if lowmem and threads > 1:
        make_error("Low-memory option is not compatible with >1 threads.")
    
    if not out.endswith((".sync", ".sync.gz")):
        make_error("Strange suffix to output file name. Exiting.")
    out_dir_file = os.path.split(out)
    if out_dir_file[0] != "" and not os.path.exists(out_dir_file[0]):
        make_error("output dir (" + out_dir_file[0] + ") not found")
    if os.path.exists(out):
        make_error(out + " already exists. Overwriting not allowed.\n")

    if len(inputs) < 2:
        make_error("At least 2 input files required")
    
    if not inputs[0].endswith((".sync.gz", ".sync")):
        make_error("file " + inputs[0] +
                   " does not have *.sync.gz or *.sync extension")
    
    if threads > mp.cpu_count():
        make_error("number of requested threads (" + str(threads) + 
                   ") exceeds number possible (" + str(mp.cpu_count()) + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Combining sync files into " + out)
    
    if lowmem:
        lowmem_read_write(inputs, out)
    else:
        highmem_read_write(inputs, out, threads)
    
    time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("... finished (" + time_str + ")")
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

