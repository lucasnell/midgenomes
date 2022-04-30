
library(tidyverse)
library(ggrepel)
library(lubridate)
library(readxl)
library(rgdal)
library(broom)
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
    filter(year(date) == 2019) %>%
    select(biotech_id, contam_pilot, n_adults, date, lake, site, location, lat, lon) %>%
    to_utm()


pools <- read_lines("~/_data/snape/space_sync_masked.names.gz")
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$biotech_id == .x][1])

freq_df <- read_tsv( "~/_data/snape/space_snape_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools), get_snape_af)) %>%
    correct_biallelic_freqs("~/_data/snape/space_sync_masked_noblanks.sync.gz")

# Kapun et al. (2021) only used SNPs > 500 bp away from each other.
freq_df <- freq_df %>%
    filter(min_dist_filter(scaff, pos, 500))


freq_mat <- freq_df %>%
    # .[,-which(colnames(freq_df) %in% samp_df$biotech_id[samp_df$contam_pilot])] %>%
    .[,-1:-3] %>%
    as.matrix() %>%
    t()


freq_pc <- prcomp(freq_mat)

freq_pc %>% summary()
freq_pc %>% str()

freq_pc_df <- freq_pc$x %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    as_tibble() %>%
    mutate(location = map_chr(sample, ~ samp_df$location[samp_df$biotech_id == .x][1]),
           lat = map_dbl(sample, ~ samp_df$lat[samp_df$biotech_id == .x][1]),
           lon = map_dbl(sample, ~ samp_df$lon[samp_df$biotech_id == .x][1]))


freq_pc_p <- freq_pc_df %>%
    ggplot(aes(PC1, PC2)) +
    geom_point() +
    geom_text_repel(aes(label = location), force = 10) +
    scale_color_viridis_d(begin = 0.2) +
    xlab(sprintf("PC1 (%.1f %% variance explained)",
                 100 * (freq_pc$sdev^2 / sum(freq_pc$sdev^2))[1])) +
    ylab(sprintf("PC2 (%.1f %% variance explained)",
                 100 * (freq_pc$sdev^2 / sum(freq_pc$sdev^2))[2])) +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_equal(xlim = c(NA, max(freq_pc_df$PC1) * 1.5)) +
    NULL
freq_pc_p

ggsave("_plots/spatial-PC.pdf", freq_pc_p, width = 6, height = 6)




