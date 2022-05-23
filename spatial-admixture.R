
library(tidyverse)
library(lubridate)
library(parallel)


options(mc.cores = max(c(1, detectCores()-2)))

source("helpers.R")

# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(ape)

# LEFT OFF: LOOK AT THESE FILES:
list.files("~/_data/admixture", "*.tsv", full.names = TRUE) %>%
    map(read_tsv, col_types = "cd") %>%
    map(~ filter(.x, llik == max(llik)))

source("plotting_funcs.R")

best1 <- c("m0_s1072639750", "m0_s1175957472", "m0_s1362936606", "m0_s1421184202", "m0_s1495316994", "m0_s1512653052", "m0_s1692003828", "m0_s1942390373", "m0_s342227430", "m0_s598692400", "m0_s784543286", "m0_s808887918", "m0_s83862104"
)

best_trees <- map(best1, ~ read.tree(sprintf("~/_data/admixture/oag_m0/%s.mltree", .x))[[2]]) %>%
    do.call(what = c)

for (i in 1:length(best_trees)) {
    cairo_pdf(sprintf("_plots/best_trees_m0/tree%i.pdf", i), width = 4, height = 3)
    plot(best_trees[[i]], no.margin = TRUE, cex = 0.5)
    dev.off()
}

# plot_resid(sprintf("~/_data/admixture/oag_m0/%s", best1[[1]]),
#            pop_order = read_lines("~/_data/snape/space_sync_masked.names.gz"))



plot_tree("~/_data/admixture/oag_m0/m0_s1072639750")
plot_tree("~/_data/admixture/oag_m0/m0_s1072639750")


cairo_pdf("_plots/admixture3.pdf", width = 8, height = 6)
plot_tree("~/_data/admixture/oag_m3/m3_s1573164154")
dev.off()

cairo_pdf("_plots/admixture3_tree.pdf", width = 6, height = 6)
read.tree("~/_data/admixture/oag_m3/m3_s1573164154.mltree")[[2]] %>%
    plot()
dev.off()






# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




library(poolfstat)
# devtools::install_github("uqrmaie1/admixtools", dependencies = TRUE)
library(admixtools)



samp_df <- read_csv("~/Box Sync/midgenomics/full-DNA-info.csv",
                    col_types = "cfcddidcccddidiDccldd") %>%
    # Convert to more accurate site names:
    mutate(site = case_when(site == "SN" ~ "Syðri Neslönd",
                            site == "KS" ~ "Kálfaströnd",
                            site == "H" ~ "Haganes",
                            site == "V" ~ "Vindbelgur",
                            site == "BR" ~ "Fellshóll",
                            TRUE ~ NA_character_),
           location = ifelse(lake == "Mývatn", site, lake)) %>%
    filter(read == 1) %>%
    filter(year(date) == 2019) %>%
    select(biotech_id, contam_pilot, n_adults, date, lake, site, location, lat, lon) %>%
    to_utm()


pools <- read_lines("~/_data/snape/space_sync_masked.names.gz") %>%
    map_chr(~ samp_df$location[samp_df$biotech_id == .x])
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$location == .x][1])

pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/space_sync_masked_noblanks.sync.gz",
                             poolsizes = n_adults * 2,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = 4)

snp_groups_mat <- min_dist_group_matrix(pool_dat@snp.info$Chromosome,
                                        pool_dat@snp.info$Position,
                                        100e3)
len <- snp_groups_mat[,2] - snp_groups_mat[,1] + 1

one_snp_chunk <- function(i) {
    chunk_i <- pooldata.subset(pool_dat,
                               snp.index = snp_groups_mat[i,1]:snp_groups_mat[i,2],
                               verbose = FALSE)
    fstats_i <- compute.fstats(chunk_i, verbose = FALSE)
    vals <- fstats_i@f2.values
    comps <- fstats_i@comparisons$F2
    out <- matrix(0, pool_dat@npools, pool_dat@npools,
                  dimnames = list(pool_dat@poolnames, pool_dat@poolnames))
    for (j in 1:nrow(comps)) {
        out[comps[j,"P1"], comps[j,"P2"]] <- vals[j,]
        out[comps[j,"P2"], comps[j,"P1"]] <- vals[j,]
    }
    return(out)
}

# Takes ~80 sec
f2_blocks <- mclapply(1:nrow(snp_groups_mat), one_snp_chunk) %>%
    simplify2array()
dimnames(f2_blocks)[[3]] <- paste0("l", len)

# Used to compare among graphs:
set.seed(7806453)
n_blocks <- dim(f2_blocks)[3]
train <- sample(1:n_blocks, round(n_blocks/2))


# ad_graphs <- find_graphs(f2_blocks, numadmix = 0, outpop = NULL,
#                           stop_gen = 100)
# winner <- ad_graphs %>% slice_min(score, with_ties = FALSE)
# winner$score[[1]]
# # [1] 151361.1


graph_options <- crossing(m = 1:10, r = 1:100)
# set.seed(375122337, "L'Ecuyer")
# graph_options$finds <- mclapply(graph_options$m, function(.m) {
#     find_graphs(f2_blocks, numadmix = .m,
#                 outpop = NULL, stop_gen = 200, verbose = FALSE)
# })


# Best score by `numadmix` (can be multiple):
best_scores <- ad_graphs_list %>%
    map(function(.x) {
        best_i <- which(.x$score == min(.x$score))
    })


oos_scores <- best_scores %>%
    imap_dfr(function(.x, .y) {
        i_graph <- ad_graphs_list[[.y]][.x[1],][["graph"]][[1]]
        res <- qpgraph(data = f2_blocks[,,train], i_graph,
                       f2_blocks_test = f2_blocks[,,-train])
        return(tibble(migr = .y, score = res$score))
    }) %>%
    arrange(score)
oos_scores


best_graph <- ad_graphs_list %>%
    .[[(oos_scores[["migr"]][[1]])]] %>%
    slice_min(score, with_ties = FALSE) %>%
    .[["graph"]] %>%
    .[[1]]

# unidentifiable_edges(best_graph)
plot_graph(best_graph, fix = FALSE)

next_best_graph <- ad_graphs_list %>%
    .[[(oos_scores[["migr"]][[2]])]] %>%
    slice_min(score, with_ties = FALSE) %>%
    .[["graph"]] %>%
    .[[1]]


fits = qpgraph_resample_multi(f2_blocks, list(best_graph, next_best_graph), nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)





# With Tjornin as outgroup:
# [1] 102516.5

plot_graph(winner$edges[[1]])

winner$edges[[1]] %>%
    edges_to_igraph() %>%
    unidentifiable_edges()
