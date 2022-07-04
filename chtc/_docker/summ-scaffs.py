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
import re
import math



def get_Ls_Ns(size_vec, thresholds):
    """
    Get Lx and Nx for all x corresponding to thresholds in `thresholds`.
    """
    
    size_vec.sort(reverse = True)
    thresholds.sort()
    
    if len(size_vec) == 1:
        lvec = [1] * len(thresholds)
        nvec = [size_vec[0]] * len(thresholds)
        return lvec, nvec
    
    lvec = [0] * len(thresholds)
    nvec = [0] * len(thresholds)
    
    if len(size_vec) == 0:
        return lvec, nvec
    
    csum = 0
    i = 0
    j = 0
    while i < len(size_vec) and j < len(thresholds):
        csum += size_vec[i]
        if csum >= thresholds[j]:
            nvec[j] = size_vec[i]
            lvec[j] = i+1
            j += 1
        i+=1
    
    return lvec, nvec





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
    N_regions = 0
    total_GC = 0
    # scaffold and contig size vectors
    ssizes = []
    csizes = []
    
    for line in fasta_file:
        if line.startswith(">"):
            ssizes.append(0)
            # If scaffolds ends with N/n for some reason, then we need this 
            # check to avoid extra zero values:
            if len(csizes) == 0 or csizes[-1] != 0:
                csizes.append(0)
        else:
            if len(ssizes) == 0:
                sys.stderr.write("FASTA doesn't start with header. Exiting.\n")
                sys.exit(1)
            # we don't want to include newline:
            line = line.rstrip()
            llen = len(line)
            ssizes[-1] += llen
            total_size += llen
            numN = line.count("n") + line.count("N")
            if numN == 0:
                csizes[-1] += llen
            elif numN == llen:
                if csizes[-1] != 0:
                    csizes.append(0)
                    N_regions += 1
            else:
                line_Nspl = [x for x in re.split("N+|n+", line) if x != ""]
                if (line[0] == "n" or line[0] == "N") and csizes[-1] != 0:
                    csizes.append(len(line_Nspl[0]))
                    N_regions += 1
                else:
                    csizes[-1] += len(line_Nspl[0])
                for i in range(1, len(line_Nspl)):
                    csizes.append(len(line_Nspl[i]))
                    N_regions += 1
                if line[-1] == "n" or line[-1] == "N":
                    csizes.append(0)
                    N_regions += 1
            total_N += numN
            total_GC += line.count("G")
            total_GC += line.count("g")
            total_GC += line.count("C")
            total_GC += line.count("c")
    
    if len(csizes) > 0 and csizes[-1] == 0:
        tmp = csizes.pop()

    fasta_file.close()
    
    contig_size = total_size - total_N
    
    sthreshs = [math.ceil(float(total_size) * x) for x in [0.5, 0.9]]
    sL, sN = get_Ls_Ns(ssizes, sthreshs)
    cthreshs = [math.ceil(float(contig_size) * x) for x in [0.5, 0.9]]
    cL, cN = get_Ls_Ns(csizes, cthreshs)
    
    print(">> Scaffold stats:")
    print("size = " + str(total_size))
    print("%GC = " + str(round(total_GC / total_size * 100, 3)))
    print("total N = " + str(total_N))
    print("N regions = " + str(N_regions))
    print("# sequences = " + str(len(ssizes)))
    print("N50 = " + str(sN[0]))
    print("L50 = " + str(sL[0]))
    print("N90 = " + str(sN[1]))
    print("L90 = " + str(sL[1]))
    print("min = " + str(ssizes[-1]))
    print("max = " + str(ssizes[0]))
    print("mean = " + str(round(sum(ssizes) / len(ssizes), 2)))
    
    print("\n\n>> Contig stats:")
    print("size = " + str(contig_size))
    print("%GC = " + str(round(total_GC / contig_size * 100, 3)))
    print("# sequences = " + str(len(csizes)))
    print("N50 = " + str(cN[0]))
    print("L50 = " + str(cL[0]))
    print("N90 = " + str(cN[1]))
    print("L90 = " + str(cL[1]))
    print("min = " + str(csizes[-1]))
    print("max = " + str(csizes[0]))
    print("mean = " + str(round(sum(csizes) / len(csizes), 2)))
    
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    sys.exit(0)

