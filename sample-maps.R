
#'
#' Maps for *Tanytarsus gracilentus* genome report.
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


save_plot <- function(n, p, w, h, ...) {
    cairo_pdf(sprintf("~/Stanford_Drive/UW/midgenomes/%s.pdf", n),
              width = w, height = h, bg = NA, ...)
    plot(p)
    dev.off()
}

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
myvatn_df <- readOGR(dsn = paste0("~/Box Sync/midgenomics/location_data/",
                                  "shapefiles/myvatn"),
                  layer = "Myvatn_WSGUTM28") |>
    tidy() |>
    rename(lon = long)
# Iceland outline is from GADM data (version 3.6; https://gadm.org/)
iceland_df <- readOGR(dsn = paste0("~/Box Sync/midgenomics/location_data/",
                                   "shapefiles/iceland"),
                 layer = "gadm36_ISL_0") |>
    tidy() |>
    rename(lon = long) |>
    to_utm() |>
    # -----------`
    # Shown below are two ways to filter this dataset, to avoid
    # plotting islands far from shore:
    # -----------`
    # 1. You can filter out islands that are very far from shore:
    # filter(!piece %in% c(90, 133, 143, 157, 215, 244, 257, 258, 260, 262))
    # 2. Filter for just the mainland:
    filter(piece == 1) |>
    identity()



# -------------`
# Sample location data
# -------------`
samp_df <- read_csv("~/Box Sync/midgenomics/full-DNA-info.csv",
                    col_types = "cfcddidcccddidiDccldd") |>
    filter(to_use == 1) |>
    filter(!is.na(site)) |>
    filter(site %in% c("KS", "BR")) |>
    group_by(site) |>
    summarize(lat = mean(lat), lon = mean(lon), .groups = "drop") |>
    add_row(site = "E5", lat = 65.597006, lon = -16.957340) |>
    # Convert to more accurate site names:
    mutate(site = case_when(site == "E5" ~ "Bolir",
                            site == "KS" ~ "Kálfaströnd",
                            site == "BR" ~ "Fellshóll",
                            TRUE ~ NA_character_)) |>
    to_utm()






# Iceland map ----

iceland_p <- iceland_df |>
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group), color = "gray20", size = 0.25,
                 fill = "gray80") +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 1.5, fill = NA, color = "black") +
    coord_equal() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL
save_plot("iceland_map", iceland_p, 4, 3)



# Myvatn map ----

myvatn_p <- myvatn_df |>
    ggplot(aes(x = lon, y = lat)) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 1.5, fill = "white", color = "black") +
    geom_polygon(aes(group = group, fill = hole), color = "gray20", size = 0.25) +
    geom_point(data = samp_df,
               size = 6, color = "gold4", fill = "gold", shape = 21) +
    geom_segment(data = tibble(lone = 406e3, lon = lone - 2e3, lat = 7280e3),
                 aes(xend = lone, yend = lat), size = 2) +
    geom_text(data = tibble(lat = 7280e3 + 300,  lon = 406e3-1e3),
              label = "2 km", vjust = 0, hjust = 0.5, size = 14 / 2.83465) +
    scale_fill_manual(values = c("lightskyblue1", "white"), guide = "none") +
    coord_equal(xlim = c(403118.1 - 100, 411932.0 + 100),
                ylim = c(7271491 - 100, 7282704 + 100)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL

save_plot("myvatn_map", myvatn_p, 4, 5)


