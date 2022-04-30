
library(tidyverse)
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)
library(poolfstat)

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


pools <- read_lines("~/_data/snape/time_sync_masked.names.gz")
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$biotech_id == .x][1])
n_pools <- length(pools)




af_df <- read_tsv("~/_data/snape/time_snape_masked_noblanks.sync.gz",
                    col_names = c("scaff", "pos", "ref", pools),
                    col_types = paste(c("cic", rep("c", n_pools)),
                                      collapse = ""))

# Kapun et al. (2021) only used SNPs > 500 bp away from each other.
af_df <- af_df %>%
    filter(min_dist_filter(scaff, pos, 500))

for (p in pools) {
    af_df[[p]] <- af_df[[p]] %>%
        str_split(":") %>%
        map_chr(~ .x[[7]]) %>%
        as.numeric()
}

# This one has strangely low covariances:
af_df <- af_df %>%
    select(-`KS-2-09-B_S48`)
filt_pools <- read_lines("~/_data/snape/time_sync_masked.names.gz") %>%
    keep(~ . != "KS-2-09-B_S48")


mean_pool <- af_df %>%
    mutate(af = rowSums(across(all_of(filt_pools))) / length(filt_pools)) %>%
    select(1:3, af)

cov_df <- tibble(biotech_id = filt_pools,
                 date = map(biotech_id,
                            ~ samp_df$date[samp_df$biotech_id == .x]) %>%
                     do.call(what = c),
                 site = map_chr(biotech_id,
                                ~ samp_df$site[samp_df$biotech_id == .x]),
                 af_cov = map_dbl(biotech_id,
                                  ~ cov(af_df[[.x]], mean_pool$af)),
                 af_cor = map_dbl(biotech_id,
                                  ~ cor(af_df[[.x]], mean_pool$af))) %>%
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
    # filter(biotech_id != "KS-2-09-B_S48") %>%
    mutate(n_adults = map_int(biotech_id, ~ samp_df$n_adults[samp_df$biotech_id == .x]),
           N = map_dbl(yr_gen, ~ tany_pop_df$N[tany_pop_df$yr_gen == .x])) %>%
    identity()

yvar <- "af_cov"
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
    ggplot(aes(af_cov, log(n_adults))) +
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


