
library(tidyverse)

#'
#' Get counts of the number of genes in each assembly's protein sets,
#' accounting for the fact that many genes contain multiple isoforms.
#'


spp_df <- tibble(spp = c("Aaegyp", "Asteph", "Bantar", "Cquinq",
                         "Cripar", "Csonor", "Ctenta", "Mdomes",
                         "Pakamu", "Ppemba", "Pstein", "Pvande",
                         "Tgraci"),
                 source = c("VectorBase", "VectorBase", "InsectBase", "VectorBase",
                            "GenBank", "here", "InsectBase", "InsectBase",
                            "InsectBase", "InsectBase", "here", "InsectBase",
                            "here"))

one_ngenes <- function(.spp, .source) {
    stopifnot(.source %in% c("VectorBase", "GenBank", "InsectBase", "here"))

    faa_lines <- read_lines(paste0("~/_data/_proteins/", .spp, "_proteins.faa.gz"))
    faa_lines <- faa_lines[grepl("^>", faa_lines)]
    head(faa_lines)


    if (.source == "VectorBase") {
        stopifnot(sum(grepl("gene=", faa_lines)) == length(faa_lines))
        ngenes <- faa_lines |>
            str_split(" \\| ") |>
            map_chr(\(x) x[grepl("^gene=", x)]) |>
            str_remove("^gene=") |>
            unique() |>
            length()
    } else if (.source == "InsectBase" || .source == "GenBank") {
        ngenes <- faa_lines |>
            str_split(" ") |>
            map_chr(\(x) x[1]) |>
            str_remove("^>") |>
            sub(pattern = "\\.[^\\.]+$", replacement = "") |>
            unique() |>
            length()
    } else {
        ngenes <- faa_lines |>
            str_remove("^>") |>
            sub(pattern = "\\.[^\\.]+$", replacement = "") |>
            unique() |>
            length()
    }

    cat(sprintf("%s\t%i\n", .spp, ngenes))

    invisible(NULL)
}



for (i in 1:nrow(spp_df)) {
    one_ngenes(spp_df$spp[[i]], spp_df$source[[i]])
}


# Run on 2023-06-06
# Aaegyp	14718
# Asteph	12705
# Bantar	10853
# Cquinq	15094
# Cripar	16522
# Csonor	18080
# Ctenta	20615
# Mdomes	14215
# Pakamu	14557
# Ppemba	17339
# Pstein	16105
# Pvande	17631
# Tgraci	15499



