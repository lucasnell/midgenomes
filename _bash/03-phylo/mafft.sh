#!/bin/bash


#'
#' Sequence alignment using MAFFT, trim using trimAl, then concatenate alignments.
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

tar -xzf ${TARGET}/prequal_filt.tar.gz -C ./


export OUT_DIR=mafft_aligns
mkdir ${OUT_DIR}


#' ============================================================================
#' ============================================================================
#'
#' Do alignments
#'
#' ============================================================================
#' ============================================================================


#' Script to run mafft to do alignment on a single gene.
#' Usage:
#'      ./one_align.sh [*.faa FILE]
cat << EOF > one_align.sh
#!/bin/bash
IN=\$1
NAME_BASE=\$(basename \$IN | sed 's/\..*//g')
OUT=./${OUT_DIR}/\${NAME_BASE}.faa
LOG=./${OUT_DIR}/\${NAME_BASE}.log
linsi --thread 1 \${IN} > \${OUT} 2> \${LOG}
EOF
chmod +x one_align.sh


#' Only show progess bar below if in an interactive job:
export inter_job="False"
if [[ $- == *i* ]]; then inter_job="True"; fi


# Takes ~5 min with 32 threads
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

if __name__ == "__main__":
    tasks = glob.glob("prequal_filt/*.faa.filtered")
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


rm -r prequal_filt one_align.sh




#' ============================================================================
#' ============================================================================
#'
#' Trim alignments
#'
#' ============================================================================
#' ============================================================================

cd ${OUT_DIR}



#' This is to create a RAxML-style partitions file to be use in
#' phylogeny construction.
#' It should look like this:
#'
#' WAG, part1 = 1-100
#' WAG, part2 = 101-384
#'
#' ... where part1/part2 is the gene name
#'
export ALIGNS_LENGTHS=trim_aligns.partition
export start=1
#'
#' Function to get unique sequence lengths in a fasta file
#' (return error if >1 value).
#'
#' Usage:
#' fa_unq_lens FILE_NAME
fa_unq_lens () {
python3 << EOF
fa_file = "$1"
lengths = []
with open(fa_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            lengths.append(0)
        else:
            lengths[-1] += len(line.rstrip())
unq_lengths = list(set(lengths))
if len(unq_lengths) > 1:
    raise AttributeError(">1 unique length in " + fa_file)
print(unq_lengths[0])
EOF
}

# Takes ~1 min
for f in *.faa; do
    g=${f%.faa}_trim.faa
    trimal -in $f -out $g -automated1
    check_exit_status "null" $?
    len=$(fa_unq_lens $g)
    end=$(( start + len - 1 ))
    echo "WAG, ${f%.faa} = ${start}-${end}" >> ../${ALIGNS_LENGTHS}
    start=$(( end + 1 ))
done

cd ..



#' ============================================================================
#' ============================================================================
#'
#' Concatenate alignments
#'
#' ============================================================================
#' ============================================================================



export CONCAT_ALIGNS=${OUT_DIR}_concat.faa

#'
#' Now run python code to concatenate all alignments together by species.
#' It assumes that all alignments contain all species.
#'
python3 << EOF
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
    genes_files = []
    # Get gene names from partition file to make sure order is the same:
    with open("trim_aligns.partition", "r") as f:
        for line in f:
            g = line.rstrip().split(" ")[1]
            gf = "mafft_aligns/" + g + "_trim.faa"
            genes_files.append(gf)
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


gzip ${CONCAT_ALIGNS}
gzip ${ALIGNS_LENGTHS}
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${CONCAT_ALIGNS}.gz ${ALIGNS_LENGTHS}.gz ${OUT_DIR}.tar.gz ${TARGET}/


rm -r ${OUT_DIR}
