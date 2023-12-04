#!/usr/bin/env python3

"""
Combine output from `summ-scaffs.py` and BUSCO into a single CSV file
to make it easier to input to a summary table.

It also makes the units for the CSV output the same as in the tables 
used to summarize the results.

This function prints to stdout.

Also only includes contig info because we didn't do any scaffolding.

Usage:
    pretty-csv.py -s [SEQ-SUMMARY FILE] -b [BUSCO-OUT FILE] [LABEL]
"""

import argparse
import sys
import os


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Make a pretty CSV")
    parser.add_argument("-s", "--summary", required = True,
                        help = "path to text file containing stdout from " +
                                "`summ-scaffs.py`")
    parser.add_argument("-b", "--busco", required = True,
                        help = "path to text file containing stdout from BUSCO")
    parser.add_argument("label", help="label to use as identifying 1st column")

    args = parser.parse_args()

    # Check files:
    if not os.path.exists(args.summary):
        sys.stderr.write(args.summary + " not found\n")
        sys.exit(1)
    if not os.path.exists(args.busco):
        sys.stderr.write(args.busco + " not found\n")
        sys.exit(1)

    # Read files:
    with open(args.summary, "rt") as summ_file:
        summ_lines = summ_file.readlines()
        summ_lines = [line.rstrip() for line in summ_lines]
    with open(args.busco, "rt") as busco_file:
        busco_lines = busco_file.readlines()
        busco_lines = [line.rstrip() for line in busco_lines]
    
    # To maintain the same order of dictionary items, even if python < 3.6 used:
    summ_cols = ("size", "contigs", "N50", "min", "max")
    busco_cols = ("BUSCO_C", "BUSCO_C-S", "BUSCO_C-D", "BUSCO_F", 
                  "BUSCO_M", "BUSCO_n")
    
    # These columns are in units of Mb so will be adjusted later:
    cols_in_mb = ("size", "N50", "max")
    # Lines containing each column's info
    summ_info = {"size": [x for x in summ_lines if x.startswith("size")][1],
                 "contigs": [x for x in summ_lines if x.startswith("# sequences")][1],
                 "N50": [x for x in summ_lines if x.startswith("N50")][1],
                 "min": [x for x in summ_lines if x.startswith("min")][1],
                 "max": [x for x in summ_lines if x.startswith("max")][1]}
    # Find just the numbers and adjust units if necessary:
    for k in summ_info.keys():
        x = summ_info[k].split(" = ")[1]
        if k in cols_in_mb:
            z = float(x) / 10**6
            x = f"{z:.6f}"
        summ_info[k] = x
    
    busco_info = {}
    busco_info["BUSCO_C"] = [s for s in busco_lines if 
                             "Complete BUSCOs" in s]
    busco_info["BUSCO_C-S"] = [s for s in busco_lines if 
                               "Complete and single-copy BUSCOs" in s]
    busco_info["BUSCO_C-D"] = [s for s in busco_lines if 
                               "Complete and duplicated BUSCOs" in s]
    busco_info["BUSCO_F"] = [s for s in busco_lines if 
                             "Fragmented BUSCOs" in s]
    busco_info["BUSCO_M"] = [s for s in busco_lines if 
                             "Missing BUSCOs" in s]
    busco_info["BUSCO_n"] = [s for s in busco_lines if 
                             "Total BUSCO groups" in s]
    # Find just the numbers:
    for k in busco_info.keys():
        x = busco_info[k][0].replace("|", "").strip().split("\t")[0]
        busco_info[k] = x
    
    HEADER =  "file,size,contigs,N50,min,max,"
    HEADER += "BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n"
    
    CSV_LINE = args.label + ","
    for k in summ_cols:
        CSV_LINE += summ_info[k] + ","
    for k in busco_cols[:-1]:
        CSV_LINE += busco_info[k] + ","
    CSV_LINE += busco_info[busco_cols[-1]]
    
    print(HEADER)
    print(CSV_LINE)
    
    sys.exit(0)


