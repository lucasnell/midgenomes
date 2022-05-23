
library(tidyverse)
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)
library(parallel)

options(mc.cores = max(1, detectCores() - 2))

source("helpers.R")

# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")


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
    filter(to_use == 1) %>%
    select(biotech_id, n_adults, date, location, lat, lon) %>%
    to_utm()


div_data <- samp_df$biotech_id %>%
    as.list() %>% set_names(samp_df$biotech_id) %>%
    map(function(x) {
        # x = samp_df$biotech_id[1]
        # rm(x, fn, xd, tlen)
        fn <- sprintf("~/_data/npstat/%s_npstat_w10000.stat.gz", x)
        xd <- read_tsv(fn, col_types = paste(c("ci", rep("d", 8)), collapse = ""),
                       progress = FALSE) %>%
            filter(length > 0) %>%
            mutate(biotech_id = x) %>%
            select(biotech_id, everything())
        return(xd)
    })




# one_div_boot <- function(xd, i) {
#     # inds <- sample(1:nrow(xd), nrow(xd), replace = TRUE, prob = xd$length)
#     inds <- sample(1:nrow(xd), nrow(xd), replace = TRUE)
#     xd_i <- xd[inds,]
#     wts_i <- xd_i$length / sum(xd_i$length)
#     out <- xd_i[1,"biotech_id"]
#     out[["rep"]] <- i
#     for (n in c("S", "Watterson", "Pi", "Tajima_D")) {
#         out[[n]] <- sum(xd_i[[n]] * wts_i, na.rm = TRUE)
#     }
#     return(out)
# }
#
# # Takes ~20 sec
# set.seed(1879340683, "L'Ecuyer")
# div_boots <- div_data %>%
#     lapply(function(xd) {
#         xdd <- mclapply(1:2000, function(i) one_div_boot(xd, i))
#         do.call(bind_rows, xdd)
#     })

#' Average across two lines of Drosophila melanogaster for the
#' "... rate of single-nucleotide mutations ... [in units of] per site per generation"
#' (Schrider et al. 2013; doi: https://doi.org/10.1534/genetics.113.151670)
mu <- 5.49e-9

div_df <- div_data %>%
    map_dfr(function(xd) {
        tlen <- sum(xd$length)
        xd %>%
            mutate(wt = length / tlen) %>%
            group_by(biotech_id) %>%
            summarize(across(S:Tajima_D, ~ sum(.x * wt, na.rm = TRUE)),
                      .groups = "drop")
    }) %>%
    left_join(samp_df, by = "biotech_id") %>%
    mutate(Ne = Watterson / (4 * mu))






div_p_df <- div_df %>%
    identity()




div_p_df %>%
    filter(yr_gen > 1977) %>%
    mutate(harm_N = map_dbl(yr_gen,
                            function(yg) {
                                n_gen <- 5
                                .N <- tany_pop_df %>%
                                    filter(yr_gen < yg, yr_gen >= yg - (n_gen * 0.5)) %>%
                                    .[["N"]]
                                return(1 / mean(1 / .N))
                            })) %>%
    arrange(yr_gen) %>%
    ggplot(aes(harm_N, Ne, color = site)) +
    geom_path() +
    geom_point() +
    scale_color_viridis_d(NULL, begin = 0.2, end = 0.7)




make_ylab <- function(yv) {
    stopifnot(yv %in% c("Watterson", "Pi", "Tajima_D", "Ne"))
    if (yv == "Watterson") z <- expression(theta[W])
    if (yv == "Pi") z <- expression(pi)
    if (yv == "Tajima_D") z <- "Tajima's D"
    if (yv == "Ne") z <- expression(N[e])
    return(z)
}
make_div_p <- function(yvar) {
    N_m <- diff(range(div_p_df[[yvar]])) / diff(range(log1p(tany_pop_df$N)))
    N_b <- - N_m * min(log1p(tany_pop_df$N)) + min(div_p_df[[yvar]])
    div_p_df %>%
        ggplot(aes(yr_gen, !!sym(yvar))) +
        geom_ribbon(data = tany_pop_df %>%
                      filter(yr_gen >= min(div_p_df$yr_gen)) %>%
                      mutate(!!sym(yvar) := log1p(N) * N_m + N_b),
                    aes(ymax = !!sym(yvar)), ymin = min(div_p_df[[yvar]]),
                  fill = "gray80") +
        geom_hline(yintercept = 0, color = "gray50", size = 1, linetype = "22") +
        # geom_line(data = tany_pop_df %>%
        #               filter(yr_gen >= min(div_p_df$yr_gen)) %>%
        #               mutate(!!sym(yvar) := log1p(N) * N_m + N_b),
        #           color = "gray60", size = 0.75) +
        geom_point(aes(color = site), size = 4, shape = 16) +
        scale_x_continuous(NULL, breaks = seq(1978, 2014, 2)) +
        scale_y_continuous(make_ylab(yvar),
                           sec.axis = sec_axis(~ (. - N_b) / N_m,
                                               "log(abundance + 1)",
                                               breaks = log1p(c(0, 100, 10e3)),
                                               labels = expression(0, 10^2, 10^4))) +
        scale_color_viridis_d(NULL, begin = 0.2, end = 0.7) +
        coord_cartesian(ylim = range(div_p_df[[yvar]])) +
        theme_classic() +
        theme(strip.text = element_blank(), legend.position = "top")
}

time_pi_p <- make_div_p("Pi")
time_theta_p <- make_div_p("Watterson")
time_D_p <- make_div_p("Tajima_D")

# ggsave("_plots/time_pi.pdf", time_pi_p, width = 6, height = 2.5, units = "in")
# ggsave("_plots/time_theta.pdf", time_theta_p, width = 6, height = 2.5, units = "in")
# ggsave("_plots/time_D.pdf", time_D_p, width = 6, height = 2.5, units = "in")


time_p <- time_pi_p +
    time_theta_p +
    time_D_p +
    plot_layout(ncol = 1, guides = "collect") &
    # theme(legend.position = "top")
    theme(legend.position = "none", axis.title.y.right = element_blank())

ggsave("_plots/time_div.pdf", time_p, width = 6, height = 6, units = "in")






myvatn_df <- readOGR(dsn = paste0("~/Box Sync/midgenomics/location_data/",
                                  "shapefiles/myvatn"),
                     layer = "Myvatn_WSGUTM28") %>%
    tidy() %>%
    rename(lon = long)
time_pal <- viridis(7, end = 0.9)[c(6, 4, 1)]
time_map_p <- myvatn_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group, fill = hole), color = "black", size = 0.1) +
    geom_point(data = distinct(samp_df, site, lat, lon),
               aes(color = site), size = 6, shape = 16) +
    geom_segment(data = tibble(lone = 405e3, lon = lone - 1e3, lat = 7280e3),
                 aes(xend = lone, yend = lat), size = 1) +
    geom_text(data = tibble(lat = 7280e3 + 200,  lon = 405e3-500),
              label = "1 km", vjust = 0, hjust = 0.5, size = 10 / 2.83465) +
    scale_fill_manual(values = c("lightblue1", "white"), guide = "none") +
    scale_color_viridis_d("Site", begin = 0.2, end = 0.7) +
    coord_equal() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          legend.justification = c(0.5, 0),
          panel.grid = element_blank()) +
    NULL

# ggsave("_plots/time_map.pdf", time_map_p, width = 2.5, height = 3, units = "in")

