#!/bin/bash


#'
#' Use codeml in PAML to estimate amino acid substitution rate from
#' concatenated alignments.
#'


. /app/.bashrc
conda activate phylo-env


export TARGET=/staging/lnell/phylo

mkdir working
cd working

export OUT_PREFIX=chir_codeml
export OUT_DIR=${OUT_PREFIX}
mkdir ${OUT_DIR}
cd ${OUT_DIR}



#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Input AA alignments and convert to PHYLIP format
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------

export CONCAT_ALIGNS=mafft_aligns_concat.faa
cp ${TARGET}/${CONCAT_ALIGNS}.gz ./ && gunzip ${CONCAT_ALIGNS}.gz

export ALIGN_PHYLIP=${CONCAT_ALIGNS/.faa/.phylip}

python3 << EOF
import glob
import sys

if __name__ == "__main__":
    aa_file = "${CONCAT_ALIGNS}"
    spp_aligns = {}
    species = []
    with open(aa_file, "r") as f:
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
    # Write all output to combined file:
    # width of output PHYLIP:
    nc = 60
    with open("${ALIGN_PHYLIP}", "w") as f:
        # b = f.write(" " + str(len(species)) + " " + str(n_aa) + " G\n")
        # b = f.write("G 1\n")
        # tmp_str = "1" * n_aa
        # for i in range(0, n_aa, nc):
        #     b = f.write(tmp_str[i:i+nc] + "\n")
        b = f.write(" " + str(len(species)) + " " + str(n_aa) + "\n")
        for spp in species:
            b = f.write(spp + "\n")
            for i in range(0, n_aa, nc):
                b = f.write(spp_aligns[spp][i:i+nc] + "\n")
    sys.exit(0)
EOF



#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#' Input ML tree and prepare for CODEML
#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------



#' This is the ML tree directly from RAxML-NG
export ML_TREE=chir_ml.tree
cp ${TARGET}/${ML_TREE} ./


#' This will be the one used in CODEML that has no branch lengths and
#' contains a few estimated point calibrations
export CM_INPUT_TREE=codeml_in_tree.nwk


R --vanilla << EOF
library(ape)
library(readr)

codeml_tr <- read.tree("${ML_TREE}")

nodeC <- getMRCA(codeml_tr, tip = c("Pstein", "Tgraci"))
nodeD <- getMRCA(codeml_tr, tip = c("Cripar", "Tgraci"))
nodeE <- getMRCA(codeml_tr, tip = c("Cmarin", "Bantar"))

codeml_tr\$edge.length <- NULL
codeml_tr\$node.label <- rep("", codeml_tr\$Nnode)
codeml_tr\$node.label[nodeC - Ntip(codeml_tr)] <- "'@2.327'"   # from Cranston et al. (2012)
codeml_tr\$node.label[nodeD - Ntip(codeml_tr)] <- "'@1.37573'" # from Cranston et al. (2012)
codeml_tr\$node.label[nodeE - Ntip(codeml_tr)] <- "'@0.76'"    # from timetree.org
write.tree(codeml_tr, "${CM_INPUT_TREE}")

EOF






#' AA substitution model:
cp /opt/conda/envs/phylo-env/dat/lg.dat ./


#' Adjusted control file from
#' https://github.com/abacus-gene/paml/blob/eb05a24a63bdeb0ce584c321057763df0f64496f/src/codeml.ctl


cat << EOF > codeml.ctl
seqfile = ${ALIGN_PHYLIP}      * sequence data filename
treefile = ${CM_INPUT_TREE}    * tree structure file name
outfile = ${OUT_PREFIX}_main.out     * main result file name

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
           * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

ndata = 1
clock = 1  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = lg.dat  * only used for aa seqs with model=empirical(_F)
           * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

model = 2
           * models for codons:
               * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
           * models for AAs or codon-translated AAs:
               * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
               * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
           * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
           * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
           * 13:3normal>0

icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0
           * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
           * AA: 0:rates, 1:separate

fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
kappa = 2  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
omega = .4 * initial or fixed omega, for codons or codon-based AAs

fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0  * different alphas for genes
ncatG = 8  * # of categories in dG of NSsites models

getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .5e-6
cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
method = 0  * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
EOF


# Takes ~25 min
codeml codeml.ctl \
    1> >(tee -a ${OUT_PREFIX}.stdout)


grep -v '^$' ${OUT_PREFIX}_main.out | grep -A1 "Substitution rate"
# Substitution rate is per time unit
# 0.090526

# # previous estimate:
# #     0.090545

gzip < ${ALIGN_PHYLIP} > ${ALIGN_PHYLIP}.gz
mv ${ALIGN_PHYLIP}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/
rm -r ${OUT_DIR}

cd ..
rm -r working
