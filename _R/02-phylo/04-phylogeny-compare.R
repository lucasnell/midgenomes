library(tidyverse)
library(ape)

# Original version without trimming or partitioning:
orig <- read.tree("_data/phylo/chir_ml.tree") |>
    ladderize()

# RAxML using merged partitions (39 in total):
merged <- read.tree(text = "(((Csonor:0.311911,(Pstein:0.183002,(((Tgraci:0.142414,((Pvande:0.026846,Ppemba:0.028594):0.097491,(Cripar:0.028996,Ctenta:0.019876):0.095380):0.031831):0.079199,(Bantar:0.113020,Cmarin:0.107146):0.041413):0.024789,Pakamu:0.118880):0.156090):0.140084):0.033043,((Cquinq:0.073456,Aaegyp:0.064305):0.056913,Asteph:0.154592):0.165120):0.200880,Mdomes:0.200880);") |>
    ladderize()

# RAxML using unmerged partitions (1,436 in total):
unmerged <- read.tree(text = "(((Asteph:0.155176,(Aaegyp:0.064525,Cquinq:0.073712):0.057113):0.165710,(Csonor:0.313034,(Pstein:0.183634,(Pakamu:0.119196,((Cmarin:0.107430,Bantar:0.113343):0.041540,(((Cripar:0.029055,Ctenta:0.019905):0.095649,(Ppemba:0.028648,Pvande:0.026897):0.097777):0.031919,Tgraci:0.142892):0.079431):0.024873):0.156621):0.140587):0.033103):0.201618,Mdomes:0.201618);") |>
    ladderize()

modeltest <- read.tree(text = "(((Asteph:0.156520,(Aaegyp:0.065208,Cquinq:0.074194):0.057944):0.168449,(Csonor:0.318813,(Pstein:0.186366,(Pakamu:0.120384,((Cmarin:0.109025,Bantar:0.114402):0.042142,(Tgraci:0.144855,((Cripar:0.029377,Ctenta:0.020046):0.097137,(Pvande:0.027186,Ppemba:0.028879):0.099397):0.032419):0.081253):0.025352):0.160455):0.144115):0.033708):0.205549,Mdomes:0.205549);") |>
    ladderize()

# merging via ModelFinder followed by ModelTest-NG, then RAxML-NG
mt_merged <- read.tree(text = "(((Csonor:0.316222,(Pstein:0.185230,(((Tgraci:0.143694,((Pvande:0.027031,Ppemba:0.028656):0.098421,(Cripar:0.029158,Ctenta:0.019919):0.096342):0.032217):0.080356,(Cmarin:0.108348,Bantar:0.113455):0.041668):0.025169,Pakamu:0.119583):0.158749):0.142484):0.034187,((Cquinq:0.073425,Aaegyp:0.064767):0.057564,Asteph:0.154777):0.166175):0.203337,Mdomes:0.203337);") |>
    ladderize()


{
    par(mfrow = c(1, 3))
    plot(orig, main = "original"); nodelabels(); tiplabels()
    plot(merged, main = "merged"); nodelabels(); tiplabels()
    plot(unmerged, main = "unmerged"); nodelabels(); tiplabels()
    # plot(modeltest, main = "modeltest"); nodelabels(); tiplabels()
}


orig_cp <- cophenetic(orig)
# merged_cp <- cophenetic(merged)[rownames(orig_cp), colnames(orig_cp)]
# unmerged_cp <- cophenetic(unmerged)[rownames(orig_cp), colnames(orig_cp)]
mt_merged_cp <- cophenetic(mt_merged)[rownames(orig_cp), colnames(orig_cp)]
# modeltest_cp <- cophenetic(modeltest)[rownames(orig_cp), colnames(orig_cp)]

mean(abs(orig_cp - unmerged_cp))
mean(abs(orig_cp - mt_merged_cp))
mean(abs(orig_cp - merged_cp))
mean(abs(unmerged_cp - merged_cp))
# mean(abs(orig_cp - modeltest_cp))
# mean(abs(unmerged_cp - modeltest_cp))


node_df <- "orig,merged,unmerged,modeltest
15,15,15,15
14,14,14,14
16,16,16,16
17,26,17,17
3,13,1,1
18,27,18,18
2,11,3,3
1,12,2,2
19,17,19,19
13,1,4,4
20,18,20,20
12,2,5,5
21,19,21,21
11,10,6,6
22,20,22,22
27,25,23,23
10,8,8,8
9,9,7,7
23,21,24,24
8,3,13,9
24,22,25,25
26,24,26,27
7,7,10,11
6,6,9,10
25,23,27,26
5,4,12,12
4,5,11,13" |>
    read_csv(col_types = "iii") |>
    arrange(orig)



orig_dn <- dist.nodes(orig)
merged_dn <- dist.nodes(merged)[node_df$merged, node_df$merged]
unmerged_dn <- dist.nodes(unmerged)[node_df$unmerged, node_df$unmerged]
# modeltest_dn <- dist.nodes(modeltest)[node_df$modeltest, node_df$modeltest]



# orig_dn[1:5, 1:5]; unmerged_dn[1:5, 1:5]; modeltest_dn[1:5,1:5]
mean(abs(orig_dn - unmerged_dn))
mean(abs(orig_dn - merged_dn))
mean(abs(unmerged_dn - merged_dn))
# mean(abs(orig_dn - modeltest_dn))
# mean(abs(unmerged_dn - modeltest_dn))



