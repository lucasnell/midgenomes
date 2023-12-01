
#'
#' Make plot of BUSCO scores for all species' transcriptomes.
#'


source("_R/00-preamble.R")


# Species whose gene predictions we made:
made_here_spp <- c("Tanytarsus gracilentus", "Parochlus steinenii",
                   "Culicoides sonorensis")

# Levels of species in order they're shown in the phylogeny:
spp_lvls <- read_csv("_data/species-names-families.csv",
                     col_types = cols(),
                     progress = FALSE)[["species"]]


#'
#' Read csv output from `cds-busco.sh` script
#'
busco_df <- "species,BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n
Aaegyp,3224,3099,125,33,28,3285
Asteph,3265,3129,136,4,16,3285
Bantar,2866,2841,25,71,348,3285
Cmarin,2970,2939,31,53,262,3285
Cquinq,3123,3079,44,15,147,3285
Cripar,3045,2771,274,43,197,3285
Csonor,2956,2370,586,50,279,3285
Ctenta,2957,2889,68,119,209,3285
Mdomes,3193,3120,73,51,41,3285
Pakamu,3109,3051,58,50,126,3285
Ppemba,3060,2970,90,68,157,3285
Pstein,3125,2915,210,26,134,3285
Pvande,3074,3009,65,77,134,3285
Tgraci,3069,2887,182,33,183,3285" |>
    read_csv(col_types = cols()) |>
    mutate(species = expand_spp(species) |>
               factor(levels = rev(spp_lvls)),
           made_here = factor(species %in% made_here_spp)) |>
    mutate(across(BUSCO_C:BUSCO_M, \(x) x / BUSCO_n * 100)) |>
    select(-BUSCO_n, -BUSCO_C) |>
    pivot_longer(`BUSCO_C-S`:BUSCO_M, names_to = "category", values_to = "perc") |>
    mutate(category = factor(category,
                             levels = paste0("BUSCO_", c("C-S","C-D","F","M")),
                             labels = c("complete+\nsingle", "complete+\nduplicated",
                                        "fragmented", "missing")))

busco_p <- busco_df |>
    ggplot(aes(perc, species)) +
    geom_col(aes(fill = made_here)) +
    xlab("Percent") +
    scale_fill_manual(NULL, values = c("gray70", "black")) +
    facet_wrap(~ category, ncol = 2, scales = "free") +
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text.y = element_text(face = "italic", color = "black"))


# save_plot("cds-busco", busco_p, w = 6.5, h = 6)
