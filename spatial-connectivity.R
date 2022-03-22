
library(tidyverse)
library(lubridate)
library(readxl)
library(rgdal)
library(broom)
library(patchwork)
library(viridis)
library(poolfstat)


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

samp_df <- read_csv("~/Box Sync/midges/full-DNA-info.csv",
                    col_types = "cfcdddiDccldd") %>%
    # Convert to more accurate site names:
    mutate(site = case_when(site == "SN" ~ "Syðri Neslönd",
                            site == "KS" ~ "Kálfaströnd",
                            site == "H" ~ "Haganes",
                            site == "V" ~ "Vindbelgur",
                            site == "BR" ~ "Fellshóll",
                            TRUE ~ NA_character_),
           location = ifelse(lake == "Mývatn", site, lake)) %>%
    filter(read == 1) %>%
    select(biotech_id, n_adults, date, lake, site, location, lat, lon) %>%
    to_utm()


pools <- c("Blik-19_S6", "Lys-19_S14", "MyBR-19_S17", "MyKS-19-B_S18",
           "MySN-19_S19")
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$biotech_id == .x][1])

pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/pilot_samples_noblanks.sync.gz",
                             poolsizes = n_adults,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = 4)

# Identity method seems to work better for dealing with low number of
# individuals in KS pool
pw_fst <- compute.pairwiseFST(pool_dat, method = "Identity",
                               nsnp.per.bjack.block = pool_dat@nsnp %/% 500)
plot(pw_fst)

# fst_mat <- pw_fst@PairwiseFSTmatrix
colnames(fst_mat) <- map_chr(colnames(fst_mat),
                             ~ samp_df$location[samp_df$biotech_id == .x][1])
rownames(fst_mat) <- map_chr(rownames(fst_mat),
                             ~ samp_df$location[samp_df$biotech_id == .x][1])

fst_mat



pw_fst@values %>%
    rownames_to_column()


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

lakes_iceland_p <- iceland_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group), color = "black", size = 0.1,
                 fill = "gray80") +
    geom_point(data = samp_df %>%
                   filter(biotech_id %in% pools, lake != "Mývatn"),
               color = "dodgerblue", size = 6) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 0.75, fill = NA, color = "black") +
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
    geom_point(data = samp_df %>%
                   filter(biotech_id %in% pools, lake == "Mývatn"),
               size = 6, color = "dodgerblue3") +
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

ggsave("~/Desktop/pilots-iceland.pdf", lakes_iceland_p, width = 6, height = 4)


ggsave("~/Desktop/pilots-myvatn.pdf", lakes_myvatn_p, width = 2.4, height = 3)

