#!/usr/bin/env python3

"""
Make the units for the CSV output from summaries of assembly steps the same
as in the tables used to summarize the results.
"""

import glob
import sys



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




if __name__ == "__main__":
    
    # Count the arguments
    cargs = len(sys.argv) - 1
    patterns = sys.argv[1:]
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
            out_csv_lines.extend(csv_lines)
    
    print("\n".join(out_csv_lines))
    print("\n")
    
    sys.exit(0)


