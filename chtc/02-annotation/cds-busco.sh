#!/bin/bash


#'
#' Run BUSCO on CDS from our assemblies and others.
#' This job was run interactively with 32 threads, 64GB RAM, and 100 GB disk.
#'


. /app/.bashrc

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=cds_busco
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export CDS_DIR=cds
cp -r ${TARGET}/${CDS_DIR}  ./
check_exit_status "copying cds directory" $?

cd ${CDS_DIR}

#' ========================================
#' If present, fix headers in Cripar because they are severely wonky.
#' ========================================
if [ -f "Cripar_cds.fasta.gz" ]; then
python3 << EOF
import re
import gzip
import sys
# Read entire FASTA file into memory:
seq_names = []
seq_seqs  = []
too_wonky = False
with gzip.open("Cripar_cds.fasta.gz","rt") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            too_wonky = "protein_id" not in line
            if not too_wonky:
                nonwonky = re.sub(".*protein_id=|\] \[.*", "", line.rstrip())
                seq_names.append(nonwonky)
                seq_seqs.append("")
        else:
            if not too_wonky:
                seq_seqs[-1] += line.rstrip()
# Now re-write it without the wonky headers:
with gzip.open("Cripar_cds.fasta.gz","wt") as fasta_file:
    for i in range(len(seq_names)):
        fasta_file.write(">" + seq_names[i] + "\n")
        fasta_file.write(seq_seqs[i] + "\n")
sys.exit(0)
EOF
fi


#' ========================================
#' Go through CDS files and only keep the longest isoform for each gene.
#' ========================================

for cds_fa in *_cds.fasta.gz; do
    species=${cds_fa%_cds.fasta.gz}
    out=${species}_cds_longest.fasta
    longest-isoforms.py ${cds_fa} > ${out}
    mv ${out} ../
    rm ${cds_fa}
    unset -v species out
done
cd ..
rm -r ${CDS_DIR}


#' ========================================
#' Python script to print pretty line of CSV for a BUSCO output:
#' ========================================
cat << EOF > pretty-busco.py
#!/usr/bin/env python3

import sys
import os

busco_name = sys.argv[1]
if len(sys.argv) > 2:
    sys.stderr.write("Only one input should be provided. Exiting.\n")
    sys.exit(1)

# Check files:
if not os.path.exists(busco_name):
    sys.stderr.write(busco_name + " not found\n")
    sys.exit(1)

# Read file:
with open(busco_name, "rt") as busco_file:
    busco_lines = busco_file.readlines()
    busco_lines = [line.rstrip() for line in busco_lines]

# To maintain the same order of dictionary items, even if python < 3.6 used:
busco_cols = ("BUSCO_C", "BUSCO_C-S", "BUSCO_C-D", "BUSCO_F",
              "BUSCO_M", "BUSCO_n")

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

CSV_LINE = ""
for k in busco_cols[:-1]:
    CSV_LINE += busco_info[k] + ","
CSV_LINE += busco_info[busco_cols[-1]]

print(CSV_LINE)

sys.exit(0)

EOF

chmod +x pretty-busco.py




#' ========================================
#' Run BUSCO on CDS and output summary to CSV
#' ========================================

conda activate busco-env

busco --download diptera_odb10 > busco-download.out
export BUSCO_DOWNLOADS=$(pwd)/busco_downloads/lineages/diptera_odb10


echo -e "species,BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n" > cds-busco.csv
for cds_fa in *_cds_longest.fasta; do

    sp=${cds_fa%_cds_longest.fasta}

    mkdir ${sp}

    BUSCO_DIR=${sp}/${sp}_busco

    busco -m transcriptome -i ${cds_fa} -o ${BUSCO_DIR} \
        --cpu ${THREADS} -l ${BUSCO_DOWNLOADS} \
    1> ${sp}/${sp}_busco.stdout \
    2> ${sp}/${sp}_busco.stderr
    check_exit_status "busco ${sp}" $?

    echo -n "${sp}," >> cds-busco.csv
    ./pretty-busco.py ${sp}/${sp}_busco.stdout >> cds-busco.csv

    unset -v BUSCO_DIR sp

done







#' ========================================
#' Manage output
#' ========================================

cd ..

cp ${OUT_DIR}/cds-busco.csv ${TARGET}/
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}
