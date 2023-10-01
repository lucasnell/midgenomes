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



absrel <- \(x) paste0("~/_data/chir_hyphy_absrel/", x)
busted <- \(x) paste0("~/_data/chir_hyphy_busted/", x)
relax <- \(x) paste0("~/_data/chir_hyphy_relax/", x)


all_absrel_data <- list.files(absrel(""), "*.json") |>
    map(\(x) { tryCatch(read_json(absrel(x)), error = function(e) NA) })
all_busted_data <- list.files(busted(""), "*.json") |>
    map(\(x) { tryCatch(read_json(busted(x)), error = function(e) NA) })
all_relax_data <- list.files(relax(""), "*.json") |>
    map(\(x) { tryCatch(read_json(relax(x)), error = function(e) NA) })
names(all_absrel_data) <- list.files(absrel(""), "*.json") |> str_remove(".json$")
names(all_busted_data) <- list.files(busted(""), "*.json") |> str_remove(".json$")
names(all_relax_data) <- list.files(relax(""), "*.json") |> str_remove(".json$")

str(all_absrel_data[[1]], max.level = 1)
str(all_busted_data[[1]], max.level = 1)
str(all_relax_data[[1]], max.level = 1)


#' The same 12 HOGs failed for all 3 No idea why.
mean(map_lgl(all_absrel_data, \(x) isTRUE(is.na(x))))
mean(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))
mean(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))
identical(names(which(map_lgl(all_absrel_data, \(x) isTRUE(is.na(x))))),
          names(which(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))))
identical(names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))),
          names(which(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))))
names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x)))))


# Remove these failed ones.
absrel_data <- all_absrel_data[map_lgl(all_absrel_data, \(x) !isTRUE(is.na(x)))]
busted_data <- all_busted_data[map_lgl(all_busted_data, \(x) !isTRUE(is.na(x)))]
relax_data <- all_relax_data[map_lgl(all_relax_data, \(x) !isTRUE(is.na(x)))]

# These objects take up a lot of memory, so remove them from environment:
rm(all_absrel_data, all_busted_data, all_relax_data); invisible(gc())


# # For aBSREL, this is how you get individual p-values:
# absrel_data[[1]][["branch attributes"]][["0"]] |>
#     map_dbl(\(x) {z <- x$`Corrected P-value`; ifelse(is.null(z), NA_real_, z)})


hyphy_df <- tibble(hog = names(absrel_data),
                   absrel_signifs = map_dbl(absrel_data, \(x) {
                       x[["test results"]][["positive test results"]]
                   }),
                   busted_pvals = map_dbl(busted_data, \(x) {
                       x[["test results"]][["p-value"]]
                   }),
                   relax_pvals = map_dbl(relax_data, \(x) {
                       x[["test results"]][["p-value"]]
                   }),
                   relax_k = map_dbl(relax_data, \(x) {
                       n <- "relaxation or intensification parameter"
                       return(x[["test results"]][[n]])
                   }),
                   busted_disc = BY_correct(busted_pvals, hog),
                   relax_disc = BY_correct(relax_pvals, hog))

table(hyphy_df$absrel_signifs)
sum(hyphy_df$busted_disc)
sum(hyphy_df$relax_disc)
sum(hyphy_df$busted_disc & hyphy_df$relax_disc)

#' From docs: "A significant result of k>1 indicates that selection strength
#' has been intensified along the test branches, and a significant result of
#' k<1 indicates that selection strength has been relaxed" along the test
#' branches.
sum(hyphy_df$relax_disc & hyphy_df$relax_k > 1)
sum(hyphy_df$relax_disc & hyphy_df$relax_k < 1)

sum(hyphy_df$busted_disc & hyphy_df$relax_disc & hyphy_df$relax_k > 1)
sum(hyphy_df$busted_disc & hyphy_df$relax_disc & hyphy_df$relax_k < 1)




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


