
library(tidyverse)
library(ggrepel)
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)

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
    filter(year(date) < 2019) %>%
    filter(to_use == 1) %>%
    select(biotech_id, n_adults, date, site)

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

count_df <- read_tsv("~/_data/snape/time_sync_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools), split_sync_strings))

freq_df <- read_tsv( "~/_data/snape/time_snape_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools), get_snape_af)) %>%
    correct_biallelic_freqs(count_df)

count_df <- count_df %>%
    filter(count_alleles(count_df %>% select(all_of(pools)) %>% as.list()) == 2)



pools2 <- read_lines("~/_data/snape/all_sync_masked.names.gz")

count_df2 <- read_tsv("~/_data/snape/all_sync_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools2),
                     col_types = paste(c("cic", rep("c", length(pools2))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools2), split_sync_strings))

freq_df2 <- read_tsv( "~/_data/snape/all_snape_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools2),
                     col_types = paste(c("cic", rep("c", length(pools2))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools2), get_snape_af)) %>%
    correct_biallelic_freqs(count_df2)

count_df2 <- count_df2 %>%
    filter(count_alleles(count_df2 %>% select(all_of(pools2)) %>% as.list()) == 2)

freq_df2 %>%
    summarize(across(all_of(pools2), ~ var(.x))) %>%
    pivot_longer(everything()) %>%
    arrange(desc(value)) %>%
    print(n = 50)



time_info_df <- samp_df %>%
    arrange(date) %>%
    mutate(gen = str_split(biotech_id, "-") %>%
               map_chr(~ .x[[2]]) %>%
               as.numeric(),
           gen = year(date) + (gen - 1) / 2,
           yr_gen = gen,
           gen = as.integer((gen - min(gen)) * 2)) %>%
    mutate(reads = lapply(biotech_id,
                            function(p) map_int(count_df[[p]], sum)))


#' NOTE:
#' From Jonas et al. (2016):
#' When the sampling variance is large compared to the drift variance ...,
#' the deviation of the outlier estimates from the true Ne is
#' particularly large.
#' For a few cases, we even observe large negative estimates.
#' Negative estimates, in general, can be interpreted as Ne being infinity,
#' that is, no evidence of genetic drift (Peel et al. 2013).
#'
calc_Ne <- function(t, pool_correct = TRUE) {
    if (t == 1) return(NA_real_)
    # rm(t, time, p0, pt, pbar, get_C, C0, Ct, F_est, Ne)
    stopifnot(length(t) == 1 && t > 1)
    t0 <- t-1
    if (time_info_df$gen[t] == time_info_df$gen[t0]) t0 <- t - 2

    time <- time_info_df$gen[t] - time_info_df$gen[t0]
    p0 <- freq_df[[time_info_df$biotech_id[t0]]]
    pt <- freq_df[[time_info_df$biotech_id[t]]]
    pbar <- (p0 + pt) / 2

    if (pool_correct) {
        get_C <- function(j) {
            Sj <- time_info_df$n_adults[j]
            Rj <- time_info_df$reads[[j]]
            Cj <- 1 / (2 * Sj) + 1 / Rj - 1 / (2 * Sj * Rj)
        }
        C0 <- get_C(t0)
        Ct <- get_C(t)
        F_est <- sum((p0 - pt)^2 - (pbar - p0 * pt) * (C0 + Ct)) /
            sum((pbar - p0 * pt) * (1 - Ct))
    } else {
        F_est <- sum((p0 - pt)^2) / sum(pbar - p0 * pt)
    }

    # Ne <- -time / (2 * log(1 - F_est))
    # This one is the "exact" version:
    Ne <- 1 / (2 * (1 - (1 - F_est)^(1 / time)))
    if (Ne < 0) Ne <- 0

    return(Ne)
}
#' For this one, it has to be t and t-1
calc_Ne_simple <- function(t) {
    if (t == 1) return(NA_real_)
    # rm(t, time, p0, pt, pbar, get_C, C0, Ct, F_est, Ne)
    stopifnot(length(t) == 1 && t > 1)

    t0 <- which(time_info_df$gen == (time_info_df$gen[t] - 1))
    if (length(t0) == 0) return(NA_real_)

    dp <- freq_df[[time_info_df$biotech_id[t]]] -
        freq_df[[time_info_df$biotech_id[t0]]]
    p <- mean(freq_df[[time_info_df$biotech_id[t]]])
    Ne <- (p * (1 - p)) / (2 * var(dp))

    return(Ne)
}
time_info_df$Ne <- map_dbl(1:nrow(time_info_df), calc_Ne)
time_info_df$Ne2 <- map_dbl(1:nrow(time_info_df), calc_Ne, pool_correct = FALSE)
time_info_df$Ne3 <- map_dbl(1:nrow(time_info_df), calc_Ne_simple)
time_info_df$pvar <- map_dbl(1:nrow(time_info_df),
                             ~ sd(freq_df[[time_info_df$biotech_id[.x]]]))

yvar <- "pvar"

p_time_info_df <- time_info_df %>%
    filter(biotech_id != "KS-2-09-B_S48") %>%
    filter(!is.na(!!sym(yvar)))

N_m <- diff(range(p_time_info_df[[yvar]], na.rm = TRUE)) / diff(range(log1p(tany_pop_df$N)))
N_b <- - N_m * min(log1p(tany_pop_df$N)) + min(p_time_info_df[[yvar]], na.rm = TRUE)

p_time_info_df%>%
    ggplot(aes(yr_gen, !!sym(yvar))) +
    geom_area(data = tany_pop_df %>%
                  filter(yr_gen >= min(p_time_info_df$yr_gen)) %>%
                  mutate(!!sym(yvar) := log1p(N) * N_m + N_b),
              fill = "gray90") +
    # geom_hline(yintercept = 50, color = "gray70", linetype = 2) +
    geom_point(aes(color = site), size = 4, shape = 16) +
    # facet_wrap(~ period, ncol = 1, scales = "free_x") +
    theme_minimal() +
    scale_x_continuous("Date", breaks = seq(1978, 2014, 2)) +
    scale_y_continuous(yvar,
                       sec.axis = sec_axis(~ (. - N_b) / N_m, "log(abundance + 1)")) +
    scale_color_viridis_d(NULL, begin = 0.2, end = 0.7) +
    coord_cartesian(ylim = range(p_time_info_df[[yvar]], na.rm = TRUE)) +
    theme(strip.text = element_blank(), legend.position = "top")








# ===========================================================================*
# ===========================================================================*

# Covariances / correlations ----

# ===========================================================================*
# ===========================================================================*


# # Kapun et al. (2021) only used SNPs > 500 bp away from each other.
# filt_freq_df <- freq_df %>%
#     filter(min_dist_filter(scaff, pos, 500))

filt_freq_df <- freq_df2 %>%
    select(1:3, all_of(pools))

# # This one has strangely low covariances:
# filt_freq_df <- filt_freq_df %>%
#     select(-`KS-2-09-B_S48`)
# filt_pools <- read_lines("~/_data/snape/time_sync_masked.names.gz") %>%
#     keep(~ . != "KS-2-09-B_S48")




# mean_pool <- filt_freq_df %>%
#     mutate(af = rowSums(across(all_of(pools))) / length(pools)) %>%
#     select(1:3, af)

mean_pool <- freq_df2 %>%
    select(1:3, `Vik-19_S41`) %>%
    rename(af = `Vik-19_S41`)

cov_df <- tibble(biotech_id = pools,
                 date = map(biotech_id,
                            ~ samp_df$date[samp_df$biotech_id == .x]) %>%
                     do.call(what = c),
                 site = map_chr(biotech_id,
                                ~ samp_df$site[samp_df$biotech_id == .x]),
                 af_cov = map_dbl(biotech_id,
                                  ~ cov(filt_freq_df[[.x]], mean_pool$af)),
                 af_cor = map_dbl(biotech_id,
                                  ~ cor(filt_freq_df[[.x]], mean_pool$af))) %>%
    mutate(gen = str_split(biotech_id, "-") %>% map_chr(~ .x[[2]]),
           gen = factor(gen, levels = paste(1:2),
                        labels = c("spring", "summer")),
           yr_gen = year(date) + ifelse(gen == "spring", 0, 0.5),
           year = factor(year(date), levels = 1977L:2015L),
           period = case_when(yr_gen < 1990 ~ 0,
                              yr_gen < 2003 ~ 1,
                              TRUE ~ 2))

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



cov_p_df <- cov_df %>%
    # filter(year(date) > 1995) %>%
    filter(biotech_id != "KS-2-09-B_S48") %>%
    mutate(n_adults = map_int(biotech_id, ~ samp_df$n_adults[samp_df$biotech_id == .x]),
           N = map_dbl(yr_gen, ~ tany_pop_df$N[tany_pop_df$yr_gen == .x])) %>%
    identity()

yvar <- "af_cor"
N_m <- diff(range(cov_p_df[[yvar]])) / diff(range(log1p(tany_pop_df$N)))
N_b <- - N_m * min(log1p(tany_pop_df$N)) + min(cov_p_df[[yvar]])


# archive_ts_p <-
cov_p_df %>%
    ggplot(aes(yr_gen, !!sym(yvar))) +
    geom_area(data = tany_pop_df %>%
                  filter(yr_gen >= min(cov_p_df$yr_gen)) %>%
                  mutate(!!sym(yvar) := log1p(N) * N_m + N_b),
              fill = "gray90") +
    # geom_hline(yintercept = 50, color = "gray70", linetype = 2) +
    geom_point(aes(color = site), size = 4, shape = 16) +
    facet_wrap(~ period, ncol = 1, scales = "free_x") +
    theme_minimal() +
    scale_x_continuous("Date", breaks = seq(1978, 2014, 2)) +
    scale_y_continuous(yvar,
                       sec.axis = sec_axis(~ (. - N_b) / N_m, "log(abundance + 1)")) +
    scale_color_viridis_d(NULL, begin = 0.2, end = 0.7) +
    coord_cartesian(ylim = range(cov_p_df[[yvar]])) +
    theme(strip.text = element_blank(), legend.position = "top")



cov_p_df %>%
    ggplot(aes(log(N), log(n_adults))) +
    # ggplot(aes(af_cov, log(N))) +
    geom_point()


lm(af_cov ~ log(N), cov_p_df) %>% {print(AIC(.)); summary(.)}
lm(af_cov ~ log(n_adults), cov_p_df) %>% {print(AIC(.)); summary(.)}
lm(af_cov ~ log(N) + log(n_adults), cov_p_df) %>% {print(AIC(.)); summary(.)}


# pool names ordered by time
ordered_pools <- cov_df %>%
    arrange(yr_gen) %>%
    .[["biotech_id"]]

pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/time_sync_masked_noblanks.sync.gz",
                             poolsizes = n_adults,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = 4)

pw_fst <- compute.pairwiseFST(pool_dat, method = "Identity",
                              nsnp.per.bjack.block = 0)
                              # nsnp.per.bjack.block = pool_dat@nsnp %/% 500)
plot(pw_fst, cex = 0.5)
heatmap(pw_fst)

fst_mat <- pw_fst@PairwiseFSTmatrix
fst_mat <- fst_mat[ordered_pools, ordered_pools]
fst_mat


