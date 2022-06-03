#!/usr/bin/env python3

"""
Make the units for the CSV output from summaries of assembly steps the same
as in the tables used to summarize the results.
"""

import glob
import sys
from pathlib import Path


def fix_one_line(in_csv_line):
    csv_list = in_csv_line.split(",")
    for i in (1, 3, 5):
        z = float(csv_list[i]) / 10**6
        csv_list[i] = f"{z:.6f}"
    csv_line = ",".join(csv_list)
    return csv_line


def fixed_csv_lines(txt):
    """ 
    Find all CSV lines in summary and 'fix' them so that the units
    (Mb in some cases) match the format used to keep track of results.
    """
    lst = txt.split("\n")
    header_idx = []
    for i, x in enumerate(lst):
        if x == HEADER:
            header_idx.append(i)
    if len(header_idx) == 0:
        return ""
    csv_idx = [x + 1 for x in header_idx]
    csv_line_list = [fix_one_line(lst[j]) for j in csv_idx]
    return csv_line_list


def rename_csv_col1(csv_lines, filename):
    """ 
    Rename first column in CSV so that it matches the file name.
    """
    new_name = Path(filename).stem
    for i in range(len(csv_lines)):
        line_i = csv_lines[i].split(",")
        line_i[0] = new_name
        csv_lines[i] = ",".join(line_i)
    return csv_lines


if __name__ == "__main__":
    
    # Count the arguments
    cargs = len(sys.argv) - 1
    patterns = sys.argv[1:]
    rename = "--rename" in patterns
    if rename:
        patterns = [x for x in patterns if x != "--rename"]
    print("\npatterns = " + " ".join(patterns) + "\n")
    if len(patterns) == 0:
        sys.exit(0)
    HEADER =  "file,size,scaffs,N50,min,max,total_N,BUSCO_C,BUSCO_C-S,BUSCO_C-D,"
    HEADER += "BUSCO_F,BUSCO_M,BUSCO_n"
    outs = []
    for pattern in patterns:
        outs.extend(glob.glob(pattern))
    out_csv_lines = []
    for filename in outs:
        with open(filename, 'r') as f:
            txt = f.read()
            csv_lines = fixed_csv_lines(txt)
            if rename:
                csv_lines = rename_csv_col1(csv_lines, filename)
            out_csv_lines.extend(csv_lines)
    
    print("\n".join(out_csv_lines))
    print("\n")
    
    sys.exit(0)


