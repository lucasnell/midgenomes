
library(tidyverse)
library(readxl)

source(".Rprofile")

samps <- read_excel("~/Box Sync/2019/midges/archived_samples.xlsx") %>%
    mutate_at(vars(year, gen, day), as.integer) %>%
    mutate(gen = factor(gen, levels = 1:2, labels = c("spring", "summer")),
           year = factor(year, levels = 1977L:2015L))

samps %>%
    group_by(year, gen, .drop = FALSE) %>%
    summarize(N = n()) %>%
    ungroup() %>%
    # Filter for generations with...
    filter(N == 0) %>%  # no sites (totally INcomplete)
    # filter(N < 2) %>%   # one site or none (partially or totally INcomplete)
    # filter(N > 0) %>%   # one site or multiple (partially or totally complete)
    # filter(N >= 2) %>%  # multiple sites (totally complete)
    # select(-N) %>%
    print(n = 100)


# How much of 78 generations is represented by each?
samps %>%
    group_by(site, year, gen, .drop = FALSE) %>%
    summarize(N = n()) %>%
    group_by(site) %>%
    summarize(N = sum(N > 0) / 78) %>%
    # summarize(N = (sum(N > 0) + (4 * 2)) / (78 + (4 * 2))) %>%
    identity()


# ---------*
# Identifying what samples to (hopefully) get from Amanda
# ---------*

# Any site is useful:
samps %>%
    group_by(year, gen, .drop = FALSE) %>%
    summarize(N = n()) %>%
    ungroup() %>%
    filter(N == 0) %>%
    arrange(year, gen) %>%
    select(gen, year) %>%
    as.data.frame() %>%
    print(row.names = FALSE)


# Any site but the one listed:
samps %>%
    group_by(year, gen, .drop = FALSE) %>%
    summarize(N = n(), site = site[1]) %>%
    ungroup() %>%
    filter(N == 1) %>%
    arrange(year, gen) %>%
    select(gen, year, site) %>%
    as.data.frame() %>%
    print(row.names = FALSE)







samps %>%
    mutate(gen = factor(2L * (as.integer(year) - 1) + as.integer(gen),
                        levels = 1L:78L)) %>%
    group_by(gen, .drop = FALSE) %>%
    summarize(SN = ifelse(sum(site == "SN") > 0, 1, 0),
              KS = ifelse(sum(site == "KS") > 0, 1, 0),
              H = ifelse(sum(site == "H") > 0, 1, 0),
              V = ifelse(sum(site == "V") > 0, 1, 0),
              sites = SN + KS + H + V) %>%
    ungroup() %>%
    # mutate(gen = paste(gen) %>% as.integer()) %>%
    gather("site", "n_samps", SN:V) %>%
    mutate(site = factor(site, levels = c("SN", "KS", "H", "V"),
                         labels = c("Syðri Neslönd", "Kálfaströnd", "Haganes", "Vindbelgur"))) %>%
    ggplot(aes(gen, n_samps, fill = site)) +
    geom_col() +
    scale_x_discrete(NULL,
                     # "Generation",
                     # breaks = seq(10, 80, 10),
                     # labels = 1977 + seq(10, 80, 10) / 2 - 1) +
                     breaks = seq(1, 79, 2),
                     labels = do.call(c, map(1977 + (seq(1, 78, 4) - 1) / 2, ~ c(.x, "")))) +
    ylab("Number of samples") +
    theme_classic() +
    theme(legend.position = "top",
          legend.direction = "vertical",
          legend.justification = c(1,0.5),
          axis.text.x = element_text(angle = 300, hjust = 0.5),
          panel.grid.major.x = element_line(size = 0.5, color = "gray80")) +
    scale_fill_brewer(NULL, palette = "Dark2")



