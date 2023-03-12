#!/bin/bash


#'
#' Ultrametric timetree based on ML tree and fossil records using mcmctree
#' in PAML.
#'
#' NOTE: PAML is mostly serial, so just `$N_MCMC_RUNS` threads is fine to
#' run the final steps in parallel.
#'

#' If you ever increase this, change `SEEDS` below.
export N_MCMC_RUNS=4

#' Ideally you provide a number of threads equal to the number of eventual
#' MCMCtree runs.
export THREADS=$(count_threads)

#' Make sure this array is at least as long as the number of replicate
#' MCMCtree runs + 1.
export SEEDS=(1920112495 2119836989 966584517 783825211 350114673 772799293
              211442557  373194819  94249041  543169895 853508965)
if (( N_MCMC_RUNS+1 > ${#SEEDS[@]} )); then
    echo "ERROR: not enough seeds provided" 1>&2
    exit 1
fi

. /app/.bashrc
conda activate phylo-env

export TARGET=/staging/lnell/phylo

mkdir working
cd working

export OUT_DIR=chir_mcmctree
mkdir ${OUT_DIR}
cd ${OUT_DIR}
export OUT_PREFIX=chir_mcmctree



export ALIGN_PHYLIP=mafft_aligns_concat.phylip
cp ${TARGET}/${ALIGN_PHYLIP}.gz ./ && gunzip ${ALIGN_PHYLIP}.gz
check_exit_status "cp, gunzip alignments" $?




#' This should be the ML tree with no branch lengths but with node calibrations
export MCMC_IN_TREE=mcmctree_in_tree.nwk
cp ${TARGET}/${MCMC_IN_TREE} ./
check_exit_status "cp calibration tree" $?




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
export BDparas_PRIORS="3.6612684 3.6612434 0.2195345"

#' Estimates for rgene_gamma priors (alpha_mu, beta_mu, alpha, respectively):
alpha=2
#' Estimate from codeml:
est_rate="0.090545"
beta=$(python -c "print('{:.5f}'.format($alpha / $est_rate))")
export rgene_gamma_PRIORS="${alpha} ${beta} 1"



#' Create control file:
cat << EOF > mcmctree.ctl
seed = ${SEEDS[0]}
seqfile = ${ALIGN_PHYLIP}
treefile = ${MCMC_IN_TREE}
mcmcfile = ${OUT_PREFIX}_mcmc.txt
outfile = ${OUT_PREFIX}.out
seqtype = 2 * 0 : nucleotides; 1: codons; 2: AAs
usedata = 3 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 1 * 1: global clock; 2: independent; and 3: correlated rates
** RootAge = '<2.954' * maximum constraint on root age
RootAge = 'U(2.954, 0.1)' * maximum constraint on root age
ndata = 1
model = 0 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0 * alpha for gamma rates at sites
ncatG = 5 * No. categories in discrete gamma
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = ${BDparas_PRIORS} * birth, death, sampling
kappa_gamma = 6 2 * gamma prior for kappa
alpha_gamma = 1 1 * gamma prior for alpha
rgene_gamma = ${rgene_gamma_PRIORS} * gammaDir prior for rate for genes
sigma2_gamma = 1 10 1 * gammaDir prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...
print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 2000
sampfreq = 5
nsample = 20000
EOF

mcmctree mcmctree.ctl \
    1> >(tee -a ${OUT_PREFIX}.stdout)

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
#' (I'm going to use the LG model instead of WAG.)


rm out.BV rst

cp /opt/conda/envs/phylo-env/dat/lg.dat ./




#' I also need to edit the "tmp0001.ctl" file. According to the docs, it
#' should look like below (with * <<<<<<< indicating changed lines).
#' Note that because `fix_alpha = 0` and `alpha > 0`, this is a discrete
#' gamma model where the estimate of alpha is estimated by codeml.

cat << EOF > tmp0001.ctl
seqfile = tmp0001.txt
treefile = tmp0001.trees
outfile = tmp0001.out
noisy = 3
seqtype = 2
model = 2 * 2: Empirical * <<<<<<<
aaRatefile = lg.dat      * <<<<<<<
fix_alpha = 0            * <<<<<<<
alpha = .5               * <<<<<<<
ncatG = 4                * <<<<<<<
Small_Diff = 0.1e-6
getSE = 2
method = 1
EOF

#'
#' From docs:
#' > This new control file will run with an empirical rate matrix
#' > (WAG [LG used here]) and with gamma rates among sites.
#' > Now you can call CODEML ... to generate the
#' > appropriate Hessian matrix using WAG+Gamma [LG+Gamma here].
#' > Rename file "rst2" as "in.BV" and now you have a nice Hessian matrix
#' > calculated using WAG+Gamma [LG+Gamma here].
#'

codeml tmp0001.ctl

# Time used: 1:08:39


mv rst2 in.BV


#' --------------------
#' Create directories for the different runs of MCMCtree, each with a
#' different starting seed.
#' --------------------

for ((i=1; i<=${N_MCMC_RUNS}; i++)); do
    MCMC_DIR=mcmc_${i}
    mkdir ${MCMC_DIR}
    sed 's/^usedata.*/usedata\ =\ 2 in.BV/g' mcmctree.ctl \
        | sed "s/^seed.*/seed\ =\ ${SEEDS[$i]}/g" \
        > ./${MCMC_DIR}/mcmctree.ctl
    cp in.BV ./${MCMC_DIR}/
    cp ${ALIGN_PHYLIP} ./${MCMC_DIR}/
    cp ${MCMC_IN_TREE} ./${MCMC_DIR}/
    unset MCMC_DIR
done


#' Simple script to run a single iteration of MCMCtree.
#' Usage:
#'  run_mcmctree.sh [RUN NUMBER]
cat << EOF > run_mcmctree.sh
#!/bin/bash
cd mcmc_\${1}
status=$?
if (( status != 0 )); then
    exit 1
fi
mcmctree mcmctree.ctl \\
    1> mcmctree.stdout \\
    2> mcmctree.stderr
status=$?
cd ..
exit $status
EOF
chmod +x run_mcmctree.sh


#'
#' Python code to run MCMCtree multiple times using multiple threads.
#'
python3 << EOF
import subprocess as sp
import multiprocessing as mp
import sys

def work(mcmc_run):
    """Defines the work unit for one run"""
    cmd = "./run_mcmctree.sh " + str(mcmc_run)
    ret = sp.run(cmd, shell = True)
    return ret

if __name__ == "__main__":
    tasks = [x+1 for x in range(${N_MCMC_RUNS})]
    n_tasks = len(tasks)
    with mp.Pool(processes=${THREADS}) as pool:
        return_codes = pool.map(work, tasks)
        return_codes = [str(x) for x in return_codes]
        print("return codes:\n" + "\n".join(return_codes))
    sys.exit(0)

EOF



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/phylo/

rm -r ${OUT_DIR}
