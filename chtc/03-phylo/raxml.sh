#!/bin/bash


#'
#' Concatenate alignments, then create phylogeny using RAxML-NG.
#'


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

export OUT_DIR=chir_raxml
export PREFIX=chir_phy
mkdir ${OUT_DIR}
cd ${OUT_DIR}



export ALIGNS_DIR=mafft_aligns
export CONCAT_ALIGNS=${ALIGNS_DIR}_concat.faa

tar -xzf ${TARGET}/${ALIGNS_DIR}.tar.gz -C ./


#'
#' Now run python code to concatenate all alignments together by species.
#' It assumes that all alignments contain all species.
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
    genes_files = glob.glob("${ALIGNS_DIR}/*.faa")
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
    with open("${CONCAT_ALIGNS}", "w") as f:
        for spp in species:
            b = f.write(">" + spp + "\n")
            for i in range(0, n_aa, nc):
                b = f.write(all_out[spp][i:i+nc] + "\n")
    sys.exit(0)

EOF

rm -r ${ALIGNS_DIR}


raxml-ng --all --msa ${CONCAT_ALIGNS} --prefix ${PREFIX} --threads ${THREADS} \
    --outgroup Anopheles_stephensi \
    --data-type AA \
    --model LG+I+G \
    --bs-trees 100 \
    --seed 453418559 \
    1> >(tee -a ${PREFIX}.stdout)


gzip ${CONCAT_ALIGNS}
mv ${CONCAT_ALIGNS}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}



