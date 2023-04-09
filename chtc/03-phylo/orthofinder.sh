#!/bin/bash


#'
#' Find orthogroups using OrthoFinder.
#' Ran this using interactive job with 48 threads, 64 GB RAM, 150 GB disk.
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_orthofinder
mkdir ${OUT_DIR}
cd ${OUT_DIR}

export PROT_FOLDER=chir_proteins

cp -r /staging/lnell/proteins ./ \
    && mv proteins ${PROT_FOLDER} \
    && cd ${PROT_FOLDER}
check_exit_status "moving, renaming protein folder" $?


# Only keep the longest isoform for each gene.
# It's assumed that the transcript names for the protein sequences
# are named as <gene name>.<transcript name> or <gene name>-<transcript name>
# and that the last '.' or '-' separates the gene name from the
# transcript name.
# It tests the first FASTA header for whether a '-' or '.' shows up last,
# and whichever is last in that header is assumed to separate gene from
# transcript names in all headers.
# If there are transcript names containing '.' or '-', then this script will
# fail in ways that may not be obvious.
#
# usage:
#     ./longest-isoforms.py in.faa > out.faa
#
cat << EOF > longest-isoform.py
#!/usr/bin/env python3
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
EOF
check_exit_status "creating longest-isoform.py" $?
chmod +x longest-isoform.py


for f in *.faa.gz; do
    g=$(echo ${f%.gz} | sed 's/_proteins//g')
    echo ${f/_proteins.faa.gz/} >> longest-isoform.log
    ./longest-isoform.py $f 1> ${g} 2>> longest-isoform.log
    status=$?
    echo -e "\n\n" >> longest-isoform.log
    if (( status != 0 )); then break; fi
    rm $f
done
check_exit_status "removing duplicate proteins" $status

cd ..


#' Create simplified time tree from MCMCTree output.
#' Use this here in OrthoFinder and save for potentially using elsewhere.
export SPECIES_TREE_MCMCTREE=chir_mcmctree.tre
export SPECIES_TREE=chir_mcmctree.nwk
cat /staging/lnell/phylo/${SPECIES_TREE_MCMCTREE} \
    | sed -e 's/\[[^][]*\]//g' \
    | grep "UTREE" \
    | sed 's/[^(]*//' \
    > ${SPECIES_TREE}
check_exit_status "create simple newick time tree" $?
cp ${SPECIES_TREE} /staging/lnell/phylo/

orthofinder -f ${PROT_FOLDER} -t ${THREADS} -a $(( THREADS / 4 )) \
    -s ${SPECIES_TREE} \
    2>&1 \
    | tee orthofinder.log
check_exit_status "run OrthoFinder" $?

# Move OrthoFinder output out of the proteins folder and rename:
cd ${PROT_FOLDER} \
    && mv OrthoFinder OrthoFinder_tmp \
    && cd OrthoFinder_tmp \
    && mv Results_* orthofinder-output \
    && mv orthofinder-output ../../ \
    && cd .. \
    && rm -r OrthoFinder_tmp \
    && cd ..
check_exit_status "move, rename OrthoFinder output" $?



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

