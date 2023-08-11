library(jsonlite)
library(tidyverse)

# LEFT OFF ----
#' How to summarize this information?
#' Are these tests redundant with each other?
#'


#' Benjamini–Yekutieli procedure for adjusting for multiple tests
#'
#' Returns a logical vector for whether each p-value represents an actual
#' discovery.
#'
BY_correct <- function(Pvals, .names, fdr = 0.05) {
    m <- length(Pvals)
    cm <- sum(1 / (1:m))
    names(Pvals) <- .names
    P <- sort(Pvals)
    no_discoveries <- TRUE
    for (k in m:1) {
        ts <- k / (m * cm) * fdr
        if (P[k] <= ts) {
            no_discoveries <- FALSE
            break
        }
    }
    if (no_discoveries) return(rep(FALSE, m))
    disc_lgl <- c(rep(TRUE, k), rep(FALSE, m-k))
    names(disc_lgl) <- names(P)
    disc_lgl <- disc_lgl[.names]
    return(disc_lgl)
}



busted <- \(x) paste0("~/_data/chir_hyphy_busted/", x)
relax <- \(x) paste0("~/_data/chir_hyphy_relax/", x)


all_busted_data <- list.files(busted(""), "*.json") |>
    map(\(x) { tryCatch(read_json(busted(x)), error = function(e) NA) })
all_relax_data <- list.files(relax(""), "*.json") |>
    map(\(x) { tryCatch(read_json(relax(x)), error = function(e) NA) })
names(all_busted_data) <- list.files(busted(""), "*.json") |> str_remove(".json$")
names(all_relax_data) <- list.files(relax(""), "*.json") |> str_remove(".json$")

str(all_busted_data[[1]], max.level = 1)
str(all_relax_data[[1]], max.level = 1)


#' The same 12 HOGs failed for both. No idea why.
mean(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))
mean(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))
identical(names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))),
          names(which(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))))
names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x)))))


# Remove these failed ones.
busted_data <- all_busted_data[map_lgl(all_busted_data, \(x) !isTRUE(is.na(x)))]
relax_data <- all_relax_data[map_lgl(all_relax_data, \(x) !isTRUE(is.na(x)))]

# These objects take up a lot of memory, so remove them from environment:
rm(all_busted_data, all_relax_data); invisible(gc())


hyphy_df <- tibble(hog = names(busted_data),
                   busted_pvals = map_dbl(busted_data, \(x) x[["test results"]][["p-value"]]),
                   relax_pvals = map_dbl(relax_data, \(x) x[["test results"]][["p-value"]]),
                   busted_disc = BY_correct(busted_pvals, hog),
                   relax_disc = BY_correct(relax_pvals, hog),
                   relax_k = map_dbl(relax_data,
                                     \(x) {
                                         z <- paste("relaxation or",
                                                    "intensification",
                                                    "parameter")
                                         return(x[["test results"]][[z]])
                                     }))

sum(hyphy_df$busted_disc)
sum(hyphy_df$relax_disc)
sum(hyphy_df$busted_disc & hyphy_df$relax_disc)



hist(hyphy_df$busted_pvals)
hist(hyphy_df$relax_pvals)


hyphy_df |>
    # filter(relax_k > 1) |>
    ggplot(aes(relax_k, - log10(relax_pvals))) +
    geom_point()


hyphy_df |>
    mutate(hog = factor(hog) |> fct_reorder(busted_pvals)) |>
    select(hog, busted_pvals, relax_pvals) |>
    ggplot(aes(hog, busted_pvals)) +
    geom_point() +
    coord_flip()

hyphy_df |>
    mutate(hog = factor(hog) |> fct_reorder(relax_pvals)) |>
    select(hog, busted_pvals, relax_pvals) |>
    ggplot(aes(hog, relax_pvals)) +
    geom_point() +
    coord_flip()


