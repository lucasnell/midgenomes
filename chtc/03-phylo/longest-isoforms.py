#!/usr/bin/env python3
"""
Only keep the longest isoform for each gene.
It's assumed that the transcript names for the protein sequences
are named as <gene name>.<transcript name> or <gene name>-<transcript name>
and that the last '.' or '-' separates the gene name from the 
transcript name.
It tests the first FASTA header for whether a '-' or '.' shows up last, 
and whichever is last in that header is assumed to separate gene from 
transcript names in all headers.
If there are transcript ids containing '.' or '-', then this script will
fail in ways that are not obvious.


usage:
    ./longest-isoforms.py in.faa > out.faa

"""


import sys
import os.path
import gzip


if len(sys.argv) > 2:
    sys.stderr.write("Only one input should be provided. Exiting.\n")
    sys.exit(1)

fasta_name = sys.argv[1]

if not os.path.exists(fasta_name):
    sys.stderr.write(fasta_name + " not found\n")
    sys.exit(1)

fasta_suffs = (".fasta.gz", ".fasta", ".fa.gz", ".fa", ".faa", ".faa.gz")

if not fasta_name.endswith(fasta_suffs):
    sys.stderr.write("Strange suffix to input FASTA file (allowed: " + 
                     ", ".join(fasta_suffs) + "). Exiting.\n")
    sys.exit(1)

if fasta_name.endswith(".gz"):
    fasta_file = gzip.open(fasta_name,"rt")
else:
    fasta_file = open(fasta_name, "r")

# Read entire FASTA into memory:
seq_names = []
seq_seqs  = []

for line in fasta_file:
    if line.startswith(">"):
        seq_names.append(line.rstrip().split(" ")[0][1:])
        seq_seqs.append("")
    else:
        seq_seqs[-1] += line.rstrip()

fasta_file.close()

# Find which type of delimeter to use:
ppos = seq_names[0].find(".")
hpos = seq_names[0].find("-")
if ppos < 0 and hpos < 0:
    sys.stderr.write("Unknown delimeter between gene and trans. names ('.' " +
                     "and '-' allowed). Exiting.\n")
    sys.exit(1)

if ppos < 0:
    deli = "-"
elif ppos > hpos:
    deli = "."
else:
    deli = "-"

genes = [x.rsplit(deli, 1)[0] for x in seq_names]
unq_genes = list(set(genes))

if len(genes) == len(unq_genes):
    sys.stderr.write("No multi-isoform genes found.")
else:
    sys.stderr.write("Removed transcripts = " + str(len(genes) - len(unq_genes)))


# Find longest isoform and write to stdout

unq_seq = ""
for i in range(len(unq_genes)):
    g = unq_genes[i]
    g_matches = [j for j, x in enumerate(genes) if x == g]
    if len(g_matches) == 0:
        sys.stderr.write("Gene '" + g + "' not found. Exiting.\n")
        sys.exit(1)
    elif len(g_matches) == 1:
        unq_seq = seq_seqs[g_matches[0]]
    else:
        g_match_lens = [len(seq_seqs[k]) for k in g_matches]
        longest = [j for j, z in enumerate(g_match_lens) if z == max(g_match_lens)][0]
        unq_seq = seq_seqs[g_matches[longest]]
    sys.stdout.write(">" + g + "\n")
    sys.stdout.write(unq_seq + "\n")


sys.exit(0)
