#!/bin/bash


#'
#' Gene family evolution analysis using CAFE.
#' Ran this using interactive job with 48 threads, 64 GB RAM, 100 GB disk.
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_cafe
mkdir ${OUT_DIR}
cd ${OUT_DIR}

#' Ultrametric species tree:
export FULL_SPECIES_TREE=time-tree.nwk
cp ${TARGET}/phylo/${FULL_SPECIES_TREE} ./
check_exit_status "moving species tree" $?


#' HOG file for N0 node from OrthoFinder:
export HOG_FILE=orthofinder-hogs-n0.tsv
cp ${TARGET}/${HOG_FILE}.gz ./  && gunzip ${HOG_FILE}.gz
check_exit_status "moving, ungzipping HOG file" $?

#' Counts file for input to CAFE
export COUNTS_FILE=orthofinder-counts-n0.tsv

#' Convert HOG file to count file for CAFE input:
R --vanilla --slave << EOF
library(readr)
count_df <- read_tsv("${HOG_FILE}", col_types = cols())
count_df <- count_df[,-which(colnames(count_df) %in%
                                   c("OG", "Gene Tree Parent Clade"))]
#' This species only had ~58% of genes match to orthogroups, so
#' I'm removing it from the analyses:
count_df <- count_df[,-which(colnames(count_df) == "Cmarin")]

spp_names <- colnames(count_df)[-1]

for (sp in spp_names) {
    z <- rep(0L, nrow(count_df))
    na_lgl <- !is.na(count_df[[sp]])
    z[na_lgl] <- sapply(strsplit(count_df[[sp]][na_lgl], ", "), length)
    count_df[[sp]] <- z
}
count_df[["Family ID"]] <- count_df[["HOG"]]
count_df[["Desc"]] <- "(null)"
count_df <- count_df[,c("Desc", "Family ID", spp_names)]
write_tsv(count_df, "${COUNTS_FILE}")
EOF
check_exit_status "convert HOG file" $?

rm ${HOG_FILE}


#' Remove Cmarin from species tree:
export SPECIES_TREE=${FULL_SPECIES_TREE%.nwk}-noCmarin.nwk
R --vanilla --slave << EOF
library(ape)
phy <- read.tree("${FULL_SPECIES_TREE}")
phy <- drop.tip(phy, "Cmarin")
write.tree(phy, "${SPECIES_TREE}")
EOF
check_exit_status "Remove Cmarin from phy" $?

rm $FULL_SPECIES_TREE


## If you later want to create lambda trees:
## export LAMBDA_TREE_CHIR=chir_sep_lambda.nwk
## export LAMBDA_TREE_CHIR_CERAT=chir_cerat_sep_lambda.nwk
## R --vanilla << EOF
## # edges for family Chironomidae
## edges2 <- c(4:11, 19:25)
## phy[["edge.length"]] <- rep(1L, length(phy[["edge.length"]]))
## edge <- phy[["edge"]]
## phy[["edge.length"]][edge[,1] %in% edges2 | edge[,2] %in% edges2] <- 2L
## write.tree(phy, "${LAMBDA_TREE_CHIR}")
##
## # edges for families Chironomidae and Ceratopogonidae
## edges2 <- c(4:12, 18:25)
## phy[["edge.length"]] <- rep(1L, length(phy[["edge.length"]]))
## edge <- phy[["edge"]]
## phy[["edge.length"]][edge[,1] %in% edges2 | edge[,2] %in% edges2] <- 2L
## write.tree(phy, "${LAMBDA_TREE_CHIR_CERAT}")
## EOF



# Run for multiple levels of k and save negative log likelihoods to choose
# which k to stick with

# This takes ~ 40 min with 48 threads

echo "k,-lnL" > nLL.csv
for i in {3..9}; do
    CAFE_OUT=cafe_k${i}
    cafe5 -t ${SPECIES_TREE} -i ${COUNTS_FILE} -c ${THREADS} -p -k $i \
        -o ${CAFE_OUT} \
        > ${CAFE_OUT}.out
    echo -n "$i," >> nLL.csv
    grep -A 1 "Inferring processes for Gamma model" ${CAFE_OUT}.out \
        | grep "^Score" \
        | sed 's/Score (-lnL): *//g' \
        >> nLL.csv
    unset -v CAFE_OUT
    echo "$i finished"
done

#$ cat nLL.csv
# k,-lnL
# 3,91696.314026833
# 4,91411.877673727
# 5,91259.849022074
# 6,91169.250635369
# 7,91112.56846517
# 8,91072.542543363
# 9,91046.328272427

#' I then looked at the following plots in R:
#' ```r
#' library(readr)
#' library(ggplot2)
#' library(dplyr)
#' x = read_csv("nLL.csv")
#' ggplot(x, aes(k, `-lnL`)) + geom_point() + geom_line()
#' ggplot(filter(x, k != 6), aes(k, `-lnL`)) + geom_point() + geom_line()
#' ```
#'

# I chose k = 8 from this bc that's when the plot started flattening out.

export k=8

mkdir other_k_vals
mv cafe_k* ./other_k_vals/

mv ./other_k_vals/cafe_k${k} cafe_k${k}_run1
mv ./other_k_vals/cafe_k${k}.out cafe_k${k}_run1.out


#' Run it with k=8 two more times because of the following from the CAFE docs:
#'
#' > Always perform multiple runs to ensure convergence, especially if
#' > multiple gamma rate categories or lambdas are used.
#'


# Takes ~15 min with 48 threads
for i in {2..3}; do
    CAFE_OUT=cafe_k${k}_run${i}
    cafe5 -t ${SPECIES_TREE} -i ${COUNTS_FILE} -c ${THREADS} -p -k $k \
        -o ${CAFE_OUT} \
        > ${CAFE_OUT}.out
    unset -v CAFE_OUT
done



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR} \
    && mv ${OUT_DIR}.tar.gz ${TARGET}/
check_exit_status "move output to target directory" $?

rm -r ${OUT_DIR}

