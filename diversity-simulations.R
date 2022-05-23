
library(tidyverse)
library(readxl)
Rcpp::sourceCpp("sim_midges.cpp")
source("helpers.R")


tany_pop_df <- read_excel(paste0("~/Box Sync/zzz-archive/2020/funding/",
                                 "ASN/Chirisl_Tany_analysis.xlsx"),
                          skip = 1,  col_names = c("year", "chir_a", "chir_b",
                                                   "chir_ab", "tany_a", "tany_b")) %>%
    select(year, starts_with("tany")) %>%
    gather("gen", "N", starts_with("tany")) %>%
    mutate(year = year + ifelse(gen == "tany_a", 0, 0.5)) %>%
    arrange(year) %>%
    select(-gen) %>%
    # To play nice with the plot of numbers of midges through time:
    rename(yr_gen = year) %>%
    mutate(period = case_when(yr_gen < 1990 ~ 0,
                              yr_gen < 2003 ~ 1,
                              TRUE ~ 2))

pools <- read_lines("~/_data/snape/time_sync_masked.names.gz")

freq_df <- read_tsv( "~/_data/snape/time_snape_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref",  pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools), get_snape_af)) %>%
    correct_biallelic_freqs("~/_data/snape/time_sync_masked_noblanks.sync.gz")
gsize <- 91518485

info77 <- read_tsv("~/_data/snape/KS-1-77_S47_snape_masked_noblanks.sync",
                   col_names = c("scaff", "pos", "ref", "sync", "snape"),
                   col_types = "ciccc") %>%
    select(-sync) %>%
    mutate(snape = str_split(snape, ":"),
           p = map_chr(snape, ~ .x[[7]]) %>% as.numeric())

# Proportion of sites segregating:
pS <- length(info77$p) / gsize


total_tany_env <- new.env(); with(total_tany_env, {
    corearea <- 0.025^2 * pi
    lakearea <- 25*10^6
    # 1000 Tanytarsus per core (500e3 per m^2) is a big generation,
    # whereas the second number is from 2012 data:
    tany_per_core <- c(peak = 1000, crash = 0.04494382)
    # Number of midges over the lake:
    lake_tany <- tany_per_core * lakearea / corearea

    #' Perhaps a bit weird, but I'll assume that 0 in the trap means ~100,000
    #' midges, and that the 2012 crash coincides with `lake_tany[["crash"]]`.
    tany_b0 <- 5e3
    tany_b1 <- lake_tany[["crash"]] / tany_pop_df$N[tany_pop_df$yr_gen == 2012]
})
tany_b0 <- total_tany_env$tany_b0
tany_b1 <- total_tany_env$tany_b1
rm(total_tany_env)

N <- tany_pop_df$N * tany_b1 + tany_b0
# 1 / mean(1 / N)

avg_het(c(info77$p, rep(0, gsize - length(info77$p))))


n_sites <- 10e3L
n_segs <- round(n_sites * pS)
P0 <- c(rep(0, n_sites - n_segs), sample(info77$p, n_segs))

#' Average across two lines of Drosophila melanogaster for the
#' "... rate of single-nucleotide mutations ... [in units of] per site per generation"
#' (Schrider et al. 2013; doi: https://doi.org/10.1534/genetics.113.151670)
mu <- 5.49e-9

het_df <- sim_het(N, P0, mu, n_reps = 12, n_threads = max(1, parallel::detectCores()-2)) %>%
    do.call(what = cbind) %>%
    as.data.frame() %>%
    as_tibble()
range(het_df[,1])

