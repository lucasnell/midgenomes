
#'
#' *Tanytarsus gracilentus* sequenced sample summary
#'

library(tidyverse)
library(lubridate)
library(readxl)
library(rgdal)
library(broom)
library(patchwork)
library(viridis)



# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")


# From decimal degrees to UTM, assuming it's Iceland and using WGS84
to_utm <- function(.df, .lat = "lat", .lon = "lon") {
    .cord.dec <- SpatialPoints(cbind(.df[[.lon]], .df[[.lat]]),
                              proj4string=CRS("+proj=longlat"))
    .cord.UTM <- spTransform(.cord.dec, CRS("+proj=utm +zone=28 ellps=WGS84"))
    .df[[.lon]] <- .cord.UTM@coords[,1]
    .df[[.lat]] <- .cord.UTM@coords[,2]
    return(.df)
}


# -------------`
# Mapping data
# -------------`
myvatn_df <- readOGR(dsn = paste0("~/Box Sync/midges/location_data/",
                                  "shapefiles/myvatn"),
                  layer = "Myvatn_WSGUTM28") %>%
    tidy() %>%
    rename(lon = long)
# Iceland outline is from GADM data (version 3.6; https://gadm.org/)
iceland_df <- readOGR(dsn = paste0("~/Box Sync/midges/location_data/",
                                   "shapefiles/iceland"),
                 layer = "gadm36_ISL_0") %>%
    tidy() %>%
    rename(lon = long) %>%
    to_utm() %>%
    # -----------`
    # Shown below are two ways to filter this dataset, to avoid
    # plotting islands far from shore:
    # -----------`
    # 1. You can filter out islands that are very far from shore:
    # filter(!piece %in% c(90, 133, 143, 157, 215, 244, 257, 258, 260, 262))
    # 2. Filter for just the mainland:
    filter(piece == 1)



# -------------`
# Sequenced sample data
# -------------`
# This is for all sequenced samples (archived and spatial)
samp_df <- read_csv("~/Box Sync/midges/full-DNA-info.csv",
                    col_types = "cfcdddiDccldd") %>%
    group_by(biotech_id, date, lake, site, n_adults, lat, lon) %>%
    summarize(gb_seq = sum(mill_seqs) * 1e6 * 150 * 1e-9,
              coverage = sum(mill_seqs) * 1e6 * 150 / 100e6,
              .groups = "drop") %>%
    # Convert to more accurate site names:
    mutate(site = case_when(site == "SN" ~ "Syðri Neslönd",
                            site == "KS" ~ "Kálfaströnd",
                            site == "H" ~ "Haganes",
                            site == "V" ~ "Vindbelgur",
                            site == "BR" ~ "Fellshóll",
                            TRUE ~ NA_character_)) %>%
    to_utm()







# ============================================================================*
# ============================================================================*

## Archived samples ----

# ============================================================================*
# ============================================================================*


archive_pal <- viridis(7, end = 0.9)[c(6, 4, 1)]

archive_df <- samp_df %>%
    filter(year(date) < 2019) %>%
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
    # # For purposes of plotting, I'll just say this is 0.5 instead of zero, so I can
    # # log transform it:
    # mutate(N = ifelse(N == 0, 0.5, N)) %>%
    # To play nice with the plot of numbers of midges through time:
    rename(yr_gen = year) %>%
    mutate(period = case_when(yr_gen < 1990 ~ 0,
                              yr_gen < 2003 ~ 1,
                              TRUE ~ 2))

N_mod <- max(log1p(tany_pop_df$N)) / 50 # max(archive_df$coverage)

archive_ts_p <- archive_df %>%
    mutate(coverage = ifelse(coverage > 60, 60, coverage)) %>%
    # This is to make sure points above coverage plotting cap don't overlap
    arrange(yr_gen, site) %>%
    split(.$yr_gen) %>%
    map_dfr(function(.d) {
        if (nrow(.d) == 1) return(.d)
        if (length(unique(.d$coverage)) > 1) return(.d)
        if (nrow(.d) > 2) stop("not programmed for > 2 sites")
        .d$coverage <- .d$coverage + 2.5 * c(-1,1)
        return(.d)
    }) %>%
    ggplot(aes(yr_gen, coverage)) +
    geom_area(data = tany_pop_df %>%
                  mutate(coverage = log1p(N) / N_mod),
              fill = "gray90") +
    geom_hline(yintercept = 50, color = "gray70", linetype = 2) +
    geom_point(data = tibble(yr_gen = map(1977:2015, ~ .x + c(0, 0.5)) %>%
                                 do.call(what = c)) %>%
                   mutate(period = case_when(yr_gen < 1990 ~ 0,
                                             yr_gen < 2003 ~ 1,
                                             TRUE ~ 2) %>%
                              factor(levels = 0:2)) %>%
                   filter(!yr_gen %in% archive_df$yr_gen),
               aes(y = 0), color = "gray60", shape = 4, size = 3) +
    geom_point(aes(color = site), size = 4, shape = 16) +
    facet_wrap(~ period, ncol = 1, scales = "free_x") +
    theme_minimal() +
    scale_x_continuous("Date", breaks = seq(1978, 2014, 2)) +
    scale_y_continuous("coverage",
                       sec.axis = sec_axis(~ . * N_mod, "log(abundance + 1)")) +
    scale_color_viridis_d("Site", begin = 0.2, end = 0.7) +
    theme(strip.text = element_blank(),
          legend.position = "none")



archive_map_p <- myvatn_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group, fill = hole), color = "black", size = 0.1) +
    geom_point(data = distinct(archive_df, site, lat, lon),
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
          legend.position = c(0.5, 1),
          legend.justification = c(0.5, 0),
          panel.grid = element_blank()) +
    NULL


archive_p <- archive_map_p +
    archive_ts_p +
    plot_annotation(tag_levels = "A") +
    plot_layout(nrow = 1, widths = c(0.5, 1))


# archive_p




# ============================================================================*
# ============================================================================*

# Other lakes ----

# ============================================================================*
# ============================================================================*


lakes_df <- samp_df %>%
    filter(year(date) == 2019) %>%
    mutate(tany_num_f = ifelse(n_adults >= 40, 1, 0) %>% factor())



lakes_iceland_p <- iceland_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group), color = "black", size = 0.1,
                 fill = "gray80") +
    geom_point(data = lakes_df %>% filter(lake != "Mývatn"),
               aes(color = tany_num_f), size = 3) +
    geom_text(data = lakes_df %>%
                  filter(n_adults < 40, lake != "Mývatn") %>%
                  arrange(lake) %>%
                  mutate(lon = lon + c(-12, 0, 0, -15, 0) * 1e3,
                         lat = lat + c(12, -10, -10, 5, -10) * 1e3),
              aes(color = tany_num_f, label = n_adults), size = 8 / 2.83465,
              vjust = 1, hjust = 0.5) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 0.75, fill = NA, color = "black") +
    scale_color_manual(values = c("firebrick1", "dodgerblue"), guide = "none") +
    coord_equal() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL


lakes_myvatn_p <- myvatn_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 0.75, fill = "white", color = "black") +
    geom_polygon(aes(group = group, fill = hole), color = "black", size = 0.1) +
    geom_point(data = lakes_df %>% filter(grepl("Mývatn", lake)),
               size = 4, color = "dodgerblue3") +
    geom_segment(data = tibble(lone = 405e3, lon = lone - 1e3, lat = 7280e3),
                 aes(xend = lone, yend = lat), size = 1) +
    geom_text(data = tibble(lat = 7280e3 + 200,  lon = 405e3-500),
              label = "1 km", vjust = 0, hjust = 0.5, size = 10 / 2.83465) +
    scale_fill_manual(values = c("lightblue1", "white"), guide = "none") +
    coord_equal(xlim = c(403118.1 - 100, 411932.0 + 100),
                ylim = c(7271491 - 100, 7282704 + 100)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL






