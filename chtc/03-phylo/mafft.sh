#!/bin/bash


#'
#' Sequence alignment using MAFFT.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

conda activate phylo-env

cp /staging/lnell/phylo/prequal_filt.tar.gz ./
tar -xzf prequal_filt.tar.gz
rm prequal_filt.tar.gz


export OUT_DIR=mafft_aligns
mkdir ${OUT_DIR}




#' Script to run mafft to do alignment on a single gene.
#' Output is `X_align.fasta` for input `X.faa`
#' Usage:
#'      ./one_align.sh [*.faa FILE]
cat << EOF > one_align.sh
#!/bin/bash
IN=\$1
NAME_BASE=\$(basename \$1 | sed 's/\..*//g')
OUT=./${OUT_DIR}/\${NAME_BASE}.fasta
LOG=./${OUT_DIR}/\${NAME_BASE}.log
linsi --thread 1 \${IN} > \${OUT} 2> \${LOG}
EOF
chmod +x one_align.sh


python3 << EOF
import glob
import subprocess as sp
import multiprocessing as mp
import sys

def work(fa_file):
    """Defines the work unit on an input file"""
    cmd = "./one_align.sh " + fa_file
    ret = sp.run(cmd, shell = True)
    return ret

tasks = glob.glob("prequal_filt/*.faa.filtered")
n_tasks = len(tasks)

with mp.Pool(processes=${THREADS}) as pool:
    for i, _ in enumerate(pool.imap_unordered(work, tasks), 1):
        sys.stderr.write('\rdone {0:%}'.format(i/n_tasks))
sys.stderr.write('\n')

sys.exit(0)

EOF


rm -r prequal_filt one_align.sh

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/



export COMBINED_FASTA=cat_aligns.fasta


#'
#' Now run python code to create a single FASTA file per species,
#' then return a vector of sorted species names so that the order is
#' consistent in the combined version.
#'
python3 << EOF
import glob
import multiprocessing as mp
import sys

def one_file(fa_file):
    """Read all species' info from a single file"""
    out = {}
    current_spp = ""
    with open(fa_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                current_spp = line.rstrip().replace(">", "").split("__")[0]
                out[current_spp] = ""
            else:
                out[current_spp] += line.rstrip()
    return out

if __name__ == "__main__":
    genes_files = glob.glob("mafft_aligns/*.fasta")
    genes_files.sort()
    species = []
    with open(genes_files[0], "r") as f:
        for line in f:
            if line.startswith(">"):
                spp = line.rstrip().replace(">", "").split("__")[0]
                species.append(spp)
    species.sort()
    # Read all files:
    all_out = {x: "" for x in species}
    for gf in genes_files:
        gout = one_file(gf)
        for spp in species:
            all_out[spp] += gout[spp]
    # These should all be the same length:
    n_aa = len(all_out[species[0]])
    if any([len(all_out[x]) != n_aa for x in species]):
        weirdos = [x for x in species if len(all_out[x]) != n_aa]
        weirdos = ", ".join(weirdos)
        raise Exception("These species had different alignment lengths: {}".format(weirdos))
    # Write all output to combined file:
    # width of output FASTA:
    nc = 60
    with open("${COMBINED_FASTA}", "w") as f:
        for spp in species:
            b = f.write(">" + spp + "\n")
            for i in range(0, n_aa, nc):
                b = f.write(all_out[spp][i:i+nc] + "\n")
    sys.exit(0)

EOF



gzip ${COMBINED_FASTA}
mv ${COMBINED_FASTA}.gz /staging/lnell/phylo/


rm -r ${OUT_DIR}
