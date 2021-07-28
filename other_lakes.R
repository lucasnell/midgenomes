library(tidyverse)
library(readxl)

source(".Rprofile")


samp_df <- read_excel("~/Box Sync/2019/midges/lakes-sampled.xlsx",
                      sheet = "samples") %>%
    mutate_if(~ is.double(.), as.integer) %>%
    mutate_if(~ is.integer(.), function(x) ifelse(is.na(x), 0L, x)) %>%
    rename_with(~ gsub("Tany", "tany", .x)) %>%
    mutate(date = as.Date(sprintf("2019-%s-%i", month, day),
                          format = "%Y-%B-%d")) %>%
    select(date, lake, sample, tany, enough_tany, confirmed, notes) %>%
    # There are two lakes named Miklavatn.
    # The one on June 8 was sampled only once and is W of the other one.
    mutate(lake = ifelse(date == as.Date("2019-06-08") & lake == "Miklavatn",
                         "Miklavatn-W", lake))



# How many lakes were sampled?
samp_df %>%
    distinct(lake) %>%
    nrow() %>%
    # Accounting for 3 sites at Myvatn
    `-`(2)


# Date ranges?
samp_df %>%
    .[["date"]] %>%
    range()



# # Which samples have > 40 Tanytarsus already?
# # >> Answer: None
# samp_df %>%
#     filter(tany > 40) %>%
#     select(date, lake, sample, tany)





# ======================================================================
# ======================================================================

# Lakes to sample ----

# ======================================================================
# ======================================================================


# Which samples were marked as having enough Tanytarsus but 40 haven't
# been yet separated?
samp_df %>%
    filter(tany < 40 & enough_tany == 1) %>%
    select(date, lake, sample, tany)



# Other lakes that seem promising
samp_df %>%
    filter(lake %in%
               # Reykjavik region:
               c("Tjornin", "Ellidavatn",
                 # Snæfellsjökull region:
                 "Lysuvatn",
                 # Between Reykjavik and Northwest region:
                 "Hredhavatn",
                 # Westfjords:
                 "Laugabolsvatn",
                 NULL
                 )) %>%
    select(date, lake, sample, tany)


# Myvatn sites
samp_df %>%
    filter(grepl("^Myvatn", lake)) %>%
    select(date, lake, sample, tany)








# ======================================================================
# ======================================================================

# Lake map ----

# ======================================================================
# ======================================================================

gps_df <- read_excel("~/Box Sync/2019/midges/lakes-sampled.xlsx", "gps_points")

get_coords <- function(.df, .which) {
    M <- as.matrix(.df[,colnames(.df)[grepl("^gps", colnames(.df))]])
    apply(M, 1,
          function(.z) {
              .z <- .z[!is.na(.z)]
              if (length(.z) == 0) return(NA_real_)
              if (length(.z) > 1) {
                  .z_dec <- .z[which(c("lat", "lon") == .which)]
                  if (.z_dec != as.integer(.z_dec)) return(.z_dec)
              }
              mean(map_dbl(.z, ~ filter(gps_df, name == .x)[[.which]]))
          })
}


map_df <- read_excel("~/Box Sync/2019/midges/lakes-sampled.xlsx",
                     "locations") %>%
    mutate(lat = get_coords(., "lat"),
           lon = get_coords(., "lon")) %>%
    mutate(date = as.Date(sprintf("2019-%s-%i", month, day),
                          format = "%Y-%B-%d")) %>%
    select(date, lake, sample, lat, lon) %>%
    # There are two lakes named Miklavatn.
    # The one on June 8 was sampled only once and is W of the other one.
    mutate(lake = ifelse(date == as.Date("2019-06-08") & lake == "Miklavatn",
                         "Miklavatn-W", lake))


ic <- map_data("world", "iceland") %>%
    rename(lon = long)


lake_map_df <- map_df %>%
    group_by(lake) %>%
    summarize(lat = median(lat, na.rm = TRUE),
              lon = median(lon, na.rm = TRUE),
              .groups = "drop")


lake_map <- ggplot(ic, aes(lon, lat)) +
    geom_polygon(fill = "white", colour = "grey50") +
    coord_quickmap() +
    geom_point(data = lake_map_df, color = "red") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank()) +
    NULL
# lake_map


# # To look at the map in Google Earth:
# lake_map_df %>%
#     rename(ID = lake) %>%
#     write_csv("~/Desktop/lakes.csv")
