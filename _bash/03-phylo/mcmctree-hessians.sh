#!/bin/bash


#'
#' Prepare for MCMCTree up to the point where we compute Hessians using CODEML.
#' Now that we have a partitioned analysis, Hessian computing will probably
#' take a while, so it makes sense to run it separately.
#'
#'

. /app/.bashrc
conda activate phylo-env


export THREADS=$(count_threads)

export TARGET=/staging/lnell/phylo

mkdir working
cd working

export OUT_DIR=chir_mcmctree_hessians
export OUT_PREFIX=chir_mcmctree_hessians


mkdir ${OUT_DIR}
cd ${OUT_DIR}


#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Partition alignments
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------


export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz

#' This was generated after merging genes into partitions and doing model
#' selection using ModelFinder.
export PARTITIONS=chir_modelfinder.part
cp ${TARGET}/${PARTITIONS} ./

export ALIGN_PHYLIP_PART=mafft_aligns_part.phylip

python3 << EOF
import glob
import sys
import os

if __name__ == "__main__":
    # First read alignment file:
    concat_file = os.environ["CONCAT_ALIGNS"]
    spp_aligns = {}
    species = []
    with open(concat_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                spp = line.rstrip().replace(">", "")
                species.append(spp)
                spp_aligns[spp] = ""
            else:
                spp_aligns[spp] += line.rstrip()
    species.sort()
    # These should all be the same length:
    n_aa = len(spp_aligns[species[0]])
    if any([len(spp_aligns[x]) != n_aa for x in species]):
        weirdos = [x for x in species if len(spp_aligns[x]) != n_aa]
        weirdos = ", ".join(weirdos)
        raise Exception("These species had different alignment lengths: {}".format(weirdos))
    #
    # Read partitions file:
    part_bounds_file = os.environ["PARTITIONS"]
    starts = []
    stops = []
    lengths = []
    with open(part_bounds_file, "r") as f:
        for line in f:
            bounds_l = line.rstrip().split(" = ")[1].split("-")
            bounds_l = [int(x) for x in bounds_l]
            # changes to indices below are to convert them from 1-based, inclusive
            # indexing to indices used in python's 'range' function that is 0-based
            # and doesn't include the stop position.
            starts.append(bounds_l[0] - 1)
            stops.append(bounds_l[1])
            lengths.append(bounds_l[1] - bounds_l[0] + 1)
    n_genes = len(starts)
    n_species = len(species)
    #
    # Write all output to a partitioned PHYLIP file:
    # width of output file:
    nc = 60
    out_file = os.environ["ALIGN_PHYLIP_PART"]
    with open(out_file, "w") as f:
        for i in range(n_genes):
            length = lengths[i]
            start = starts[i]
            stop = stops[i]
            b = f.write(" " + str(n_species) + " " + str(length) + "\n")
            for spp in species:
                b = f.write(spp + "\n")
                nl = 0
                for x in range(start, stop):
                    b = f.write(spp_aligns[spp][x])
                    nl += 1
                    if nl >= nc:
                        f.write("\n")
                        nl = 0
                if nl < nc:
                    b = f.write("\n")
            if i % 100 == 0:
                print("%6.2f%%" % (i / n_genes * 100))
    sys.exit(0)
EOF


export N_GENES=$(grep "^ " $ALIGN_PHYLIP_PART | wc -l)

# Note: partitions will be used later
rm $CONCAT_ALIGNS
unset -v CONCAT_ALIGNS




#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Add calibrations to ML tree
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------


#' This is the ML tree directly from RAxML-NG
export ML_TREE=chir_ml.tree
cp ${TARGET}/${ML_TREE} ./


#' This will be the one used in MCMCTree that has no branch lengths and
#' with node calibrations
export MCMC_INPUT_TREE=mcmctree_in_tree.nwk


R --vanilla << EOF
library(ape)
library(readr)

mcmc_tr <- read.tree("${ML_TREE}")

nodeB <- getMRCA(mcmc_tr, tip = c("Csonor", "Tgraci"))
nodeC <- getMRCA(mcmc_tr, tip = c("Pstein", "Tgraci"))
nodeD <- getMRCA(mcmc_tr, tip = c("Cripar", "Tgraci"))
nodeE <- getMRCA(mcmc_tr, tip = c("Cmarin", "Bantar"))

mcmc_tr\$edge.length <- NULL
mcmc_tr\$node.label <- rep("", mcmc_tr\$Nnode)

#' The extra complexity below is for two reasons:
#'
#' 1. \`ape\` removes parentheses and commas when writing to a newick file, so
#'    I'm assigning each character a "safe" alternative, then replacing them
#'    in the output string from \`ape\` before writing to file.
#' 2. I need to add the number of species and trees at the top of the input
#'    file to mcmctree.
#'

mcmc_tr\$node.label[nodeB - Ntip(mcmc_tr)] <- "'L__2.420-0.10-0.2-0.1--'"  # B
mcmc_tr\$node.label[nodeC - Ntip(mcmc_tr)] <- "'L__2.013-0.16-0.5--'"      # C
mcmc_tr\$node.label[nodeD - Ntip(mcmc_tr)] <- "'L__0.935-0.47-0.5--'"      # D
mcmc_tr\$node.label[nodeE - Ntip(mcmc_tr)] <- "'L__0.339-1.24-1.0--'"      # E

#' Now write to newick file for use in mcmctree
mcmc_tr_file <- "${MCMC_INPUT_TREE}"
write_file(paste0(Ntip(mcmc_tr), " 1\n"), mcmc_tr_file)
mcmc_tr_str <- gsub("__", "(", write.tree(mcmc_tr))
mcmc_tr_str <- gsub("--", ")", mcmc_tr_str)
mcmc_tr_str <- gsub("-", ", ", mcmc_tr_str)
write_file(paste0(mcmc_tr_str, "\n"), mcmc_tr_file, append = TRUE)

EOF


rm ${ML_TREE}
unset -v ML_TREE


#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Setup MCMCTree parameters and run it
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

#' Note: The PAML documentation says: "The time unit should be chosen such
#' that the node ages are roughly in the range 0.01-10."
#' So I'll be using 100 Ma as the unit.

#' NOTE: See "Tutorial 4" (p 12) at http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf
#' for how to analyze protein data


#' ------------------------------------
#' Prior specifications (others below aren't used):
#' ------------------------------------
#'
#' BDparas: get relative timetree by inputting protein data and ML tree from
#'          RAxML-NG to RelTime-ML in MEGA without any calibrations
#'          (see https://www.megasoftware.net/web_help_11/Part_I_Getting_Started/A_Walk_Through_MEGA/Constructing_a_Timetree_(ML).htm).
#'          Then use this tree as input to ddBD tree prior method
#'          (see https://github.com/cathyqqtao/ddBD-tree-prior
#'          and https://doi.org/10.1093/bioinformatics/btab307).
#'          Use `root.time = 2.2`, which is the approximate estimate for this
#'          node on timetree.org in units of 100 Ma.
#' rgene_gamma: First get an estimate for the overall rate by running codeml
#'              on the ML tree using only a few point calibration estimates
#'              (e.g., change "B(0.6, 0.8)" to "0.7") and clock = 1.
#'              From the mcmctree documentation:
#'              >
#'              > The gamma distribution has mean alpha/beta and variance
#'              > alpha/beta.
#'              > The first parameter (alpha) controls the shape of the
#'              > distribution.
#'              > Values of alpha = 1 or = 2 lead to fairly diffuse priors.
#'              > It is advisable to set alpha to one of those two values, and
#'              > then fix beta so that the mean rate is reasonable.
#'              > We use a value of alpha_D = 1 for the Dirichlet parameter
#'              > [this is the 3rd option for this prior].
#'              > This is the default value, that produces a reasonable
#'              > partitioning of the rate across loci.
#'              >
#'              Thus, for an estimated rate X, the prior below should be
#'              2 2/X 1
#'

#' Estimates for BDparas priors (lambda [birth rate], mu [death rate], and
#' rho [sampling fraction], respectively):
export BDparas_PRIORS="1.0683690 1.0683750 0.7739223"


#' Estimates for rgene_gamma priors (alpha_mu, beta_mu, alpha, respectively):
alpha=2
#' Estimate from codeml:
est_rate="0.087005"
beta=$(python -c "print('{:.5f}'.format($alpha / $est_rate))")
export rgene_gamma_PRIORS="${alpha} ${beta} 1"
unset -v est_rate beta

#' We're using the correlated rates model (clock = 3).
#' For the priors for the gamma distribution used to generate sigma^2 values
#' for among-branch variability, we're using the methods from "The Impact of the
#' Representation of Fossil Calibrations on Bayesian Estimation of Species
#' Divergence Times".
#' In this paper, they set the alpha parameter to m / s and
#' the beta parameter to m / s^2, where m and s are the mean and stdev of
#' the gamma distribution, respectively.
#' Because they set m = s = 1 / rt_mean (rt_mean is the mean of the distribution
#' used to generate root age), beta simplifies to rt_mean.
#' For rt_mean, we'll use the middle of the uniform portion of the distribution
#' used to generate root ages (i.e., (max + min) / 2).
#'
export sigma2_gamma_PRIORS="1 $(python -c "print('{:.4f}'.format((2.385 + 2.954) / 2))") 1"




mkdir hessians

#' Create control file:
cat << EOF > hessians/mcmctree.ctl
seed = 2074462960
seqfile = $(pwd)/${ALIGN_PHYLIP_PART}
treefile = $(pwd)/${MCMC_INPUT_TREE}
mcmcfile = ${OUT_PREFIX}_mcmc.txt
outfile = ${OUT_PREFIX}.out
seqtype = 2 * 0 : nucleotides; 1: codons; 2: AAs
usedata = 3 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 3 * 1: global clock; 2: independent; and 3: correlated rates
RootAge = 'B(2.385, 2.954, 0.025, 0.1)' * minimum and maximum constraints on root age
ndata = ${N_GENES}
model = 0 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0 * alpha for gamma rates at sites
ncatG = 5 * No. categories in discrete gamma
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = ${BDparas_PRIORS} * birth, death, sampling
kappa_gamma = 6 2 * gamma prior for kappa
alpha_gamma = 1 1 * gamma prior for alpha
rgene_gamma = ${rgene_gamma_PRIORS} * gammaDir prior for rate for genes
sigma2_gamma = ${sigma2_gamma_PRIORS} * gammaDir prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...
print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 2000
sampfreq = 5
nsample = 20000
EOF

cd hessians

# Takes ~ 5 min
mcmctree mcmctree.ctl \
    1> ${OUT_PREFIX}.stdout



#' From documentation:
#'
#' > MCMCTree will generate the "tmp0001.ctl", "tmp0001.trees", "tmp0001.out"
#' > and "tmp0001.txt" files and will call CODEML to generate the Hessian
#' > matrix for the protein data. However, MCMCTree will use the simplest
#' > protein model (Poisson and no gamma rates) which is not very useful for
#' > real data analysis.
#' > Delete the "out.BV" and "rst" files that were generated.
#' > Copy file "wag.dat" from the "dat" directory into [your working directory].
#'
#' (I'm going to use the separate models by partition instead of WAG.)


rm out.BV rst



#' Associative array for model files for all unique best-fitting substitution
#' models (via ModelTest-NG) across all partitions.
#'
#' I made this list by manually checking the unique models from ModelTest-NG
#' output, then choosing the file inside `/opt/conda/envs/phylo-env/dat/`
#' that fit closest.
#' Because LG was the best model in 71% of models, we replaced the models
#' that don't have an equivalent in MCMCTree (PMB and VT) with LG.
#'
declare -A mod_files=( ["DAYHOFF"]="dayhoff.dat" ["DCMUT"]="dayhoff-dcmut.dat" \
                       ["JTT"]="jones.dat" ["JTT-DCMUT"]="jones-dcmut.dat" \
                       ["LG"]="lg.dat" ["MTART"]="mtArt.dat" \
                       ["MTZOA"]="MtZoa.dat" ["PMB"]="lg.dat" \
                       ["VT"]="lg.dat" ["WAG"]="wag.dat" )



#' I also need to edit the "tmp0001.ctl" file. According to the docs, it
#' should look like below (with * <<<<<<< indicating changed lines).
#' Note that because `fix_alpha = 0` and `alpha > 0`, this is a discrete
#' gamma model where the estimate of alpha is estimated by codeml.

## cat << EOF > tmp0001.ctl
## seqfile = tmp0001.txt
## treefile = tmp0001.trees
## outfile = tmp0001.out
## noisy = 3
## seqtype = 2
## Small_Diff = 0.1e-6
## getSE = 2
## method = 1
## model = 2 * 2: Empirical * <<<<<<<
## aaRatefile = lg.dat      * <<<<<<<
## fix_alpha = 0            * <<<<<<<
## alpha = .5               * <<<<<<<
## ncatG = 4                * <<<<<<<
## EOF
##


#' Make these changes:
for ((i=1; i<=N_GENES; i++)); do
    # get model and model dat file for ith partition:
    mod=$(sed -n "$i"p ../${PARTITIONS} | sed 's/+.*//g') \
        && dat="${mod_files[$mod]}"
    check_exit_status "null" $?
    # move files into separate folder:
    new_dir=part_$(printf "%04d\n" $i) \
        && mkdir ${new_dir} \
        && mv tmp$(printf "%04d\n" $i)* ./${new_dir}/ \
        && cd ${new_dir} \
        && cp /opt/conda/envs/phylo-env/dat/${dat} ./
    check_exit_status "null" $?
    # edit control file:
    ctl=tmp$(printf "%04d\n" $i).ctl \
        && sed -i "/^model/d; /^aaRatefile/d" ${ctl} \
        && echo "model = 2 * 2: Empirical" >> ${ctl} \
        && echo "aaRatefile = ${dat}" >> ${ctl} \
        && echo "fix_alpha = 0" >> ${ctl} \
        && echo "alpha = .5" >> ${ctl} \
        && echo "ncatG = 4" >> ${ctl}
    check_exit_status "null" $?
    cd ..
done


#'
#' From docs:
#' > This new control file will run with an empirical rate matrix
#' > (WAG [LG used here]) and with gamma rates among sites.
#' > Now you can call CODEML ... to generate the
#' > appropriate Hessian matrix using WAG+Gamma [LG+Gamma here].
#' > Rename file "rst2" as "in.BV" and now you have a nice Hessian matrix
#' > calculated using WAG+Gamma [LG+Gamma here].
#'

# Without partitioning, this is what was run:
# # Takes ~ 100 min
# codeml tmp0001.ctl
# mv rst2 in.BV



#' Simple script to run a single iteration of CODEML.
#' Usage:
#'  run_codeml.sh [PARTITION NUMBER]
cat << EOF > run_codeml.sh
#!/bin/bash
# add zero-padding:
pad_i=\$(printf "%04d\n" \${1})
cd part_\${pad_i}
status=\$?
if (( status != 0 )); then
    exit 1
fi
codeml tmp\${pad_i}.ctl \\
    1> codeml_\${pad_i}.stdout \\
    2> codeml_\${pad_i}.stderr
status=\$?
cd ..
exit \$status
EOF
chmod +x run_codeml.sh


#'
#' Python code to run CODEML multiple times using multiple threads.
#'
#' Takes ~3 min with 48 threads
#'
python3 << EOF
import subprocess as sp
import multiprocessing as mp
import sys
import os

def work(partition):
    """Defines the work unit for one run"""
    cmd = "./run_codeml.sh " + str(partition)
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    # one partition per gene:
    n_parts = int(os.environ["N_GENES"])
    tasks = [x+1 for x in range(n_parts)]
    n_tasks = len(tasks)
    n_threads = int(os.environ["THREADS"])
    with mp.Pool(processes=n_threads) as pool:
        return_codes = pool.map(work, tasks)
        return_codes = [str(x) for x in return_codes]
        print("\n-----------------------\n")
        print("return codes:")
        print("\n-----------------------\n")
        print("\n".join(return_codes))
    sys.exit(0)

EOF




cd ../../
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

cd ..
rm -r working
