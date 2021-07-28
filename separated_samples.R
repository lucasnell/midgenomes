
library(tidyverse)
library(readxl)

source(".Rprofile")


samples <- read_excel("~/Box Sync/2020/midges/tany_from_archive.xlsx",
                        sheet = "data-entry") %>%
    mutate(across(ends_with(c("_year", "_month", "_day")), as.integer)) %>%
    mutate(gen = factor(gen, levels = 1:2, labels = c("spring", "summer")),
           year = factor(col_year, levels = 1977L:2015L),
           date = as.Date(paste(col_year, col_month, col_day, sep = "-")),
           N = as.integer(num_males)) %>%
    select(site, year, gen, date, N)


samples %>%
    filter(N > 0) %>%
    arrange(year, gen, site, date) %>%
    group_by(site, year, gen) %>%
    mutate(n1 = sum(N), n2 = max(N)) %>%
    ungroup() %>%
    # distinct(site, year, gen, n1, n2) %>%
    filter(n2 < 40, n1 >= 10, n1 != n2)


    # summarize(N = sum(N), .groups = "drop") %>%
    group_by(site, year, gen) %>%
    # filter(N == max(N)) %>%
    summarize(N = sum(N), .groups = "drop") %>%
    ungroup()





separated <- samples %>%
    filter(N >= 20)


# Generations NOT found in `separated`
crossing(year = factor(1977L:2015L, levels = 1977L:2015L),
         gen = factor(1:2, levels = 1:2, labels = c("spring", "summer"))) %>%
    split(interaction(.$year, .$gen)) %>%
    map_dfr(function(.x) {
        .z <- filter(.x, !year %in% separated$year[separated$gen == gen])
        if (nrow(.z) > 0) {
            .zz <- samples %>%
                filter(year == .x$year, gen == .x$gen)
            if (nrow(.zz) > 0) {
                mutate(.z, N = sum(.zz[["N"]]))
            } else mutate(.z, N = 0)
        } else mutate(.z, N = integer(0))
    }) %>%
    arrange(year, gen) %>%
    identity()
    # write_csv("~/Desktop/samples.csv")



midge_df <- read_excel("~/Box Sync/2020/funding/ASN/Chirisl_Tany_analysis.xlsx",
                       skip = 1,  col_names = c("year", "chir_a", "chir_b", "chir_ab",
                                                "tany_a", "tany_b")) %>%
    select(year, starts_with("tany")) %>%
    gather("gen", "N", starts_with("tany")) %>%
    mutate(year = year + ifelse(gen == "tany_a", 0, 0.5)) %>%
    arrange(year) %>%
    select(-gen) %>%
    # For purposes of plotting, I'll just say this is 0.5 instead of zero, so I can
    # log transform it:
    mutate(N = ifelse(N == 0, 0.5, N)) %>%
    # To play nice with the plot of numbers of midges through time:
    rename(date = year) %>%
    mutate(N = (log(N) - min(log(N))) * 40 / diff(range(log(N))))


# Numbers of males through time:
needed_midges <- crossing(year = factor(1977L:2015L, levels = 1977L:2015L),
                          gen = factor(1:2, levels = 1:2,
                                       labels = c("spring", "summer")),
                          site = sort(unique(samples$site[samples$N > 0]))) %>%
    split(interaction(.$year, .$gen, .$site)) %>%
    map_dfr(function(.x) {
        .z <- samples %>%
            filter(year == .x$year, gen == .x$gen, site == .x$site)
        .N <- case_when(sum(.z$N) > 0 ~ sum(.z$N),
                        .x$site == "SN" ~ 0L,
                        TRUE ~ NA_integer_)
        mutate(.x, N = .N)
    }) %>%
    mutate(date = as.numeric(paste(year)) + ifelse(gen == "spring", 0, 0.5),
           N = case_when(N > 40 ~ 40.0,
                         TRUE ~ N * 1.0),
           Z = factor(as.integer(N > 0), levels = 0:1)) %>%
    split(interaction(.$year, .$gen, drop = TRUE)) %>%
    map_dfr(function(.x) {
        stopifnot(nrow(.x) > 0)
        logl <- duplicated(.x$N) & !is.na(.x$N) & .x$N > 0
        if (any(logl)) {
            .x$N[logl] <- .x$N[logl] + (1:sum(logl))
        }
        return(.x)
    })

sampled_p <- needed_midges %>%
    mutate(period = case_when(date < 1990 ~ 0,
                              date < 2003 ~ 1,
                              TRUE ~ 2) %>%
               factor(levels = 0:2)) %>%
    ggplot(aes(date, N)) +
    geom_vline(data = tibble(date = map(1977:2015, ~ .x + c(0, 0.5)) %>%
                                 do.call(what = c)) %>%
                   mutate(period = case_when(date < 1990 ~ 0,
                                             date < 2003 ~ 1,
                                             TRUE ~ 2) %>%
                              factor(levels = 0:2)),
               aes(xintercept = date), color = "gray70") +
    geom_area(data = midge_df %>%
                  mutate(period = case_when(date < 1990 ~ 0,
                                            date < 2003 ~ 1,
                                            TRUE ~ 2) %>%
                             factor(levels = 0:2)),
              fill = "red", alpha = 0.25) +
    geom_point(aes(shape = Z, color = site), na.rm = TRUE) +
    # geom_line(aes(color = site), na.rm = TRUE) +
    facet_wrap(~ period, nrow = 3, scales = "free_x") +
    scale_shape_manual(values = c(1, 19), guide = FALSE) +
    scale_color_manual(values = c("dodgerblue", "firebrick", "gray30")) +
    scale_x_continuous(breaks = seq(1978, 2014, 2)) +
    ylab("Number of individuals") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_blank())
# sampled_p

# ggsave("~/Desktop/midges_sampled.pdf", sampled_p, width = 6, height = 6)



needed_midges %>%
    mutate(date = factor(date)) %>%
    ggplot(aes(site, N)) +
    # geom_area(data = midge_df, fill = "red", alpha = 0.25) +
    geom_point(aes(shape = Z, color = site), na.rm = TRUE) +
    facet_wrap(~ date) +
    scale_shape_manual(values = c(1, 19), guide = FALSE) +
    scale_color_manual(values = c("dodgerblue", "firebrick", "gray30")) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.5),
          strip.background = element_blank())



library(lubridate)

other_sites <- read_excel("~/Box Sync/2019/midges/archived_samples.xlsx") %>%
    mutate_at(vars(year, gen, day), as.integer) %>%
    mutate(gen = factor(gen, levels = 1:2, labels = c("spring", "summer")),
           year = factor(year, levels = 1977L:2015L),
           date = paste(year, str_sub(month, 1, 3), day, sep = "-") %>%
               as_date(format = "%Y-%b-%d")) %>%
    filter(site != "SN") %>%
    select(year, gen, site, date)


needed_midges %>%
    filter(N < 30) %>%
    split(interaction(.$year, .$gen, drop = TRUE)) %>%
    map_dfr(function(.x) {
        z <- other_sites %>%
            filter(year == .x[["year"]], gen == .x[["gen"]])
        if (nrow(z) == 0) return(NULL)
        return(z)
    }) %>%
    arrange(year, gen, site, date) %>%
    mutate(year = paste(year), gen = paste(gen),
           year = ifelse(year == lag(year,1,"-1"), "", year),
           gen = ifelse(gen == lag(gen,1,"-1"), "", gen),
           date = paste(date)) %>%
    # identity()
    write_csv("~/Desktop/other_sites.csv")




