library(tidyverse)

#' THIS IS NOW DEFUNCT, EVERYTHING USEFUL THAT'S DONE HERE IS DONE IN `make-genome-stats.R`

repeats_df <- read_csv("_data/repeats-summary.csv", col_types = cols())

# repeats_df$length |> log10() |> hist()

#' To convert back to full names:
spp_name_map <- list("Aaegyp" = "Aedes aegypti",
                     "Asteph" = "Anopheles stephensi",
                     "Bantar" = "Belgica antarctica",
                     "Cripar" = "Chironomus riparius",
                     "Ctenta" = "Chironomus tentans",
                     "Cmarin" = "Clunio marinus",
                     "Cquinq" = "Culex quinquefasciatus",
                     "Csonor" = "Culicoides sonorensis",
                     "Mdomes" = "Musca domestica",
                     "Pstein" = "Parochlus steinenii",
                     "Ppemba" = "Polypedilum pembai",
                     "Pvande" = "Polypedilum vanderplanki",
                     "Pakamu" = "Propsilocerus akamusi",
                     "Tgraci" = "Tanytarsus gracilentus")

gsize_map <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    mutate(species = map_chr(species, \(x) spp_name_map[[x]])) |>
    select(species, gsize) |>
    (\(x) {
        z <- as.list(x$gsize)
        names(z) <- x$species
        return(z)
    })()

repeats_df |>
    mutate(species = map_chr(species, \(x) spp_name_map[[x]]),
           perc = 100 * length / map_dbl(species, \(x) gsize_map[[x]])) |>
    write_csv("_data/repeats-summary-2.csv")
