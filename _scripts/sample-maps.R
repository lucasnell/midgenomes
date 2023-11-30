
#'
#' Maps of *Tanytarsus gracilentus* sample locations.
#'

source("_scripts/00-preamble.R")

library(sf)



#' ===========================================================================
#' ===========================================================================
# Myvatn map ----
#' ===========================================================================
#' ===========================================================================

myvatn_map <- st_read("_data/spatial/Myvatn_WSGUTM28.geojson")

samp_df <- read_csv("_data/spatial/Tgraci-sample-locations.csv", col_types = "cdd")

myvatn_bounds <- st_bbox(myvatn_map) |> as.list()

myvatn_scale_df <- tibble(x = myvatn_bounds$xmin + 500,
                          x_end = x + 2e3,
                          y = myvatn_bounds$ymax - 3e3,
                          y_end = y,
                          x_mid = (x + x_end) / 2)

myvatn_p <- myvatn_map |>
    ggplot() +
    geom_sf(color = "gray20", linewidth = 0.25, fill = "lightskyblue1") +
    geom_point(data = samp_df, aes(x, y), size = 6, color = "gold4",
               fill = "gold", shape = 21) +
    geom_segment(data = myvatn_scale_df,
                 aes(x, y, xend = x_end, yend = y_end), linewidth = 2) +
    geom_text(data = myvatn_scale_df,
              aes(x_mid, y + 300), label = "2 km", vjust = 0, hjust = 0.5,
              size = 14 / 2.83465) +
    geom_rect(data = as_tibble(myvatn_bounds),
              aes(xmin = xmin - 100, xmax = xmax + 100,
                  ymin = ymin - 100, ymax = ymax + 100),
              linewidth = 1.5, fill = NA, color = "black") +
    theme_void() +
    coord_sf(ylim = unlist(myvatn_bounds[c("ymin", "ymax")]) + c(-200, 200),
             xlim = unlist(myvatn_bounds[c("xmin", "xmax")]) + c(-200, 200)) +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA))

# save_plot("myvatn_map", myvatn_p, 4, 5, .png = FALSE)





#' ===========================================================================
#' ===========================================================================
# Iceland map ----
#' ===========================================================================
#' ===========================================================================


# Iceland outline is from GADM data (version 4.1; https://gadm.org/)
iceland_map <- st_read("_data/spatial/gadm41_ISL_0.geojson") |>
    st_transform(st_crs("+proj=utm +zone=28 ellps=WGS84"))
geom_nrows <- map_int(iceland_map$geometry[[1]], \(x) nrow(x[[1]]))

st_geometry(iceland_map) <- iceland_map$geometry[[1]][[which(geom_nrows == max(geom_nrows))]] |>
    st_polygon() |>
    st_sfc(crs = st_crs(iceland_map))


iceland_p <- iceland_map |>
    ggplot() +
    geom_sf(color = "gray20", linewidth = 0.25, fill = "gray80") +
    geom_rect(data = as_tibble(myvatn_bounds),
              aes(xmin = xmin - 100, xmax = xmax + 100,
                  ymin = ymin - 100, ymax = ymax + 100),
              linewidth = 1.5, fill = NA, color = "black") +
    theme_void() +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent"))

# save_plot("iceland_map", iceland_p, 4, 3, .png = FALSE)



