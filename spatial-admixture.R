
library(tidyverse)
library(lubridate)
library(parallel)
library(poolfstat)
# devtools::install_github("uqrmaie1/admixtools", dependencies = TRUE)
library(admixtools)

options(mc.cores = max(c(1, detectCores()-2)))

source("helpers.R")





# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")


samp_df <- read_csv("~/Box Sync/midges/full-DNA-info.csv",
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


pools <- read_lines("~/_data/snape/space_sync_masked.names.gz")
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$biotech_id == .x][1])

pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/space_sync_masked_noblanks.sync.gz",
                             poolsizes = n_adults * 2,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = 4)

chroms <- unique(pool_dat@snp.info$Chromosome)
pool_dat@snp.info$Position[which(pool_dat@snp.info$Chromosome == chroms[2])[1]]
# pool_dat@snp.info$Chromosome
# pool_dat@snp.info$Position


z <- admixtools::example_f2_blocks
str(z)

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

ad_graphs <- find_graphs(f2_blocks, numadmix = 0, outpop = NULL,
                          stop_gen = 100)
winner <- ad_graphs %>% slice_min(score, with_ties = FALSE)
winner$score[[1]]
# [1] 151361.1



# With Tjornin as outgroup:
# [1] 102516.5

plot_graph(winner$edges[[1]])

winner$edges[[1]] %>%
    edges_to_igraph() %>%
    unidentifiable_edges()
