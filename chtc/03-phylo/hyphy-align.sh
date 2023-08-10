#!/bin/bash


#'
#' Do alignments for later selection tests using HyPhy.
#' Ran this using interactive job with 16 threads, 32 GB RAM, 100 GB disk.
#'


. /app/.bashrc

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_hyphy_inputs
mkdir ${OUT_DIR}
cd ${OUT_DIR}

export CDS_DIR=cds
cp -r ${TARGET}/${CDS_DIR}  ./
check_exit_status "copying cds directory" $?

cd ${CDS_DIR}

#' If present, fix headers in Cripar because they are severely wonky.
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



#' Go through CDS files and only keep the longest isoform for each gene.
#'
#' Before processing files, remove this species's file because it had
#' poor matching on OrthoFinder:
if [ -f "Cmarin_cds.fasta.gz" ]; then rm Cmarin_cds.fasta.gz; fi
# Now do the actual removing of isoforms:
for cds_fa in *_cds.fasta.gz; do
    species=${cds_fa%_cds.fasta.gz}
    out=${species}_cds_longest.fasta
    longest-isoforms.py ${cds_fa} > ${out}
    rm ${cds_fa}
    unset -v species out
done
cd ..




#' This contains all the genes for each HOG by species. It should be sorted
#' by species, gene, and HOG (in that order).
export HOG_GENES_CSV=hyphy-hog-genes.csv
cp ${TARGET}/${HOG_GENES_CSV}.gz  ./ \
    && gunzip ${HOG_GENES_CSV}.gz
check_exit_status "copying, gunzipping HOG genes CSV file" $?


# ------------
# Use `$HOG_GENES_CSV` to organize into CDS fasta files by HOG:
# -----------
export HOG_CDS_DIR=hog_${CDS_DIR}
mkdir ${HOG_CDS_DIR}
cd ${HOG_CDS_DIR}


python3 << EOF

import csv
import sys
import gzip

class HogGeneCsvInfo:
    """Extract info from CSV file of species, genes, and HOGs"""
    def __init__(self, csv_filename):
        self.species = []
        self.genes = {}
        self.hogs = {}
        with open(csv_filename, "rt") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            header = next(csv_reader, None)
            if header[0] != "species" or header[1] != "gene" or header[2] != "hog":
                err_msg = "header of CSV should be 'species', 'gene', 'hog' (in "
                err_msg += "that order). Yours is " + ", ".join(header)
                sys.stderr.write(err_msg)
                sys.exit(1)
            for row in csv_reader:
                if len(self.species) == 0 or row[0] != self.species[-1]:
                    self.species.append(row[0])
                    self.genes[self.species[-1]] = [row[1]]
                    self.hogs[self.species[-1]] = [row[2]]
                else:
                    self.genes[self.species[-1]].append(row[1])
                    self.hogs[self.species[-1]].append(row[2])
        return None
    def nspp(self):
        return len(self.species)
    def ngenes(self, sp):
        return len(self.genes[sp])
    def nhogs(self, sp):
        return len(self.hogs[sp])

class FastaInfo:
    """Extract an entire FASTA file into memory"""
    def __init__(self, filename):
        self.items = {}
        latest = ""
        if filename.endswith(".gz"):
            fasta_file = gzip.open(filename,"rt")
        else:
            fasta_file = open(filename, "r")
        for line in fasta_file:
            if line.startswith(">"):
                latest = line.rstrip().split(" ")[0][1:]
                self.items[latest] = ""
            else:
                self.items[latest] += line.rstrip()
        fasta_file.close()
        return None
    def __getitem__(self, key):
        try:
            out = self.items[key]
        except KeyError:
            out = ""
        return out

if __name__ == "__main__":
    csv_filename = "../${HOG_GENES_CSV}"
    csv_info = HogGeneCsvInfo(csv_filename)
    cds_dir = "../${CDS_DIR}/"
    for sp in csv_info.species:
        cds_filename = cds_dir + sp + "_cds_longest.fasta"
        cds_fa = FastaInfo(cds_filename)
        for j in range(csv_info.ngenes(sp)):
            g = csv_info.genes[sp][j]
            h = csv_info.hogs[sp][j]
            hog_filename = h + ".fasta"
            with open(hog_filename, "a") as hog_file:
                hog_file.write(">" + sp + "\n")
                fa_string = cds_fa[g]
                if fa_string == "":
                    err_msg = "could not find key " + g + " for species " + sp
                    sys.stderr.write(err_msg + "\n")
                    sys.exit(1)
                hog_file.write(fa_string)
                hog_file.write("\n")
    sys.exit(0)

EOF



cd ..

rm -r ${CDS_DIR}


export ALIGNS_DIR=hyphy_aligns
mkdir ${ALIGNS_DIR}


cd ${HOG_CDS_DIR}






#' Script to run hyphy's codon-msa and mafft to do alignment on a single HOG.
#' Usage:
#'      ./one_align.sh [*_cds.fasta FILE]
cat << EOF > one_align.sh
#!/bin/bash

. /app/.bashrc

IN=\$1
NAME_BASE=\${IN%.fasta}

mkdir \${NAME_BASE}
cd \${NAME_BASE}

cp ../\${IN} ./

log_file=\${NAME_BASE}.log
prot_fas=\${NAME_BASE}_protein.fas
nuc_fas=\${NAME_BASE}_nuc.fas
prot_msa=\${prot_fas%.fas}.msa
nuc_msa=\${NAME_BASE}.msa

hyphy /opt/codon-msa/pre-msa.bf --input \${IN} > \${log_file} 2>&1

mv \${IN}_protein.fas \${prot_fas} \\
    && mv \${IN}_nuc.fas \${nuc_fas}

conda activate phylo-env
linsi --thread 1 \${prot_fas} > \${prot_msa} 2>> \${log_file}
conda deactivate

hyphy /opt/codon-msa/post-msa.bf --protein-msa \${prot_msa} \\
    --nucleotide-sequences \${nuc_fas} --output \${nuc_msa} --compress No \\
     >> \${log_file} 2>&1

rm \${IN}

cd ..
mv \${NAME_BASE} ../${ALIGNS_DIR}/

exit 0
EOF
chmod +x one_align.sh


#' ---------------------------------------
#' Now use python's multiprocessing package to run hyphy's codon-msa and mafft
#' across all HOGs using multiple threads.
#' ---------------------------------------

#' Only show progess bar below if in an interactive job:
export inter_job="False"
if [[ $- == *i* ]]; then inter_job="True"; fi

python3 << EOF
import glob
import subprocess as sp
import multiprocessing as mp
import sys

def work(cds_file):
    """Defines the work unit on an input file"""
    cmd = "./one_align.sh " + cds_file
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    tasks = glob.glob("*.fasta")
    n_tasks = len(tasks)
    show_progress = $inter_job
    with mp.Pool(processes=${THREADS}) as pool:
        for i, _ in enumerate(pool.imap_unordered(work, tasks), 1):
            if show_progress:
                sys.stderr.write('\rdone {0:%}'.format(i/n_tasks))
    if show_progress:
        sys.stderr.write('\n')
    sys.exit(0)

EOF

rm one_align.sh

cd ../..

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}
