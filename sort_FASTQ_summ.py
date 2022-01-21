#!/usr/bin/env python3
"""
Sort Summary Information for a FASTQ file.

This script sorts a summary file so that the reads are in the same order as 
they are in a FASTQ file.

usage:
    ./sort_FASTQ_summ.py <summary file> <FASTQ file>

"""


import sys
import os.path
import numpy as np
import pandas as pd
import gzip





if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("You should only provide summary and FASTQ files as arguments")
        sys.exit(1)
    summ_fn = str(sys.argv[1])
    fastq_fn = str(sys.argv[2])
    if not os.path.exists(summ_fn):
        print(summ_fn + " not found")
        sys.exit(1)
    if not os.path.exists(fastq_fn):
        print(fastq_fn + " not found")
        sys.exit(1)

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
        with gzip.open(fastq_fn,"rt") as in_fastq:
            for i, line in enumerate(in_fastq):
                if i % 4 == 0:
                    if j >= len(fq_read_names):
                        print("summary file is missing read names")
                        sys.exit(1)
                     # removes all after first space
                    just_name = line.partition(" ")[0]
                    #  `[1:]` removed starting @
                    fq_read_names[j] = just_name[1:]
                    j+=1
    else:
        with open(fastq_fn, "r") as in_fastq:
            for i, line in enumerate(in_fastq):
                if i % 4 == 0:
                    if j >= len(fq_read_names):
                        print("summary file is missing read names")
                        sys.exit(1)
                    just_name = line.partition(" ")[0]
                    fq_read_names[j] = just_name[1:]
                    j+=1
    
    sort_summ_df = summ_df.reindex(fq_read_names)
    # check for whether all read names match:
    if sort_summ_df["batch_id"].isnull().values.any():
        print(sort_summ_df["batch_id"].isnull().sum(), " read names don't match")
        sys.exit(1)

    if all(sort_summ_df.index == fq_read_names):
        print("sorted read names should match those from FASTQ file!")
        sys.exit(1)
    
    # Turning index back to column:
    sort_summ_df.reset_index(inplace = True)
    # Reordering to original column order:
    sort_summ_df = sort_summ_df[col_names]
    
    # Change just the file name (keep directory the same) for output
    out_dir_file = os.path.split(summ_fn)
    out_name = os.path.join(out_dir_file[0], "sorted_" + out_dir_file[1])
    
    print("writing to " + out_name)
    
    summ_df.to_csv(out_name, sep="\t", index = False)
    
    sys.exit(0)


