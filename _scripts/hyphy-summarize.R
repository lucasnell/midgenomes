
source("_scripts/00-preamble.R")

library(jsonlite)
library(gt)
library(knitr)

# To have tibbles display more digits:
options(pillar.sigfig = 4)



#' Benjamini–Yekutieli procedure for adjusting for multiple tests
#'
#' Returns a logical vector for whether each p-value represents an actual
#' discovery.
#'
BY_correct <- function(Pvals, .names, fdr = 0.1) {
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


busted <- \(x) paste0(dirs$hyphy_busted, "/", x)
relax <- \(x) paste0(dirs$hyphy_relax, "/", x)


all_busted_data <- list.files(busted(""), "*.json") |>
    map(\(x) { tryCatch(read_json(busted(x)), error = function(e) NA) })
all_relax_data <- list.files(relax(""), "*.json") |>
    map(\(x) { tryCatch(read_json(relax(x)), error = function(e) NA) })
names(all_busted_data) <- list.files(busted(""), "*.json") |> str_remove(".json$")
names(all_relax_data) <- list.files(relax(""), "*.json") |> str_remove(".json$")


#' The same 12 HOGs failed for both.
mean(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))
mean(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))
identical(names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x))))),
          names(which(map_lgl(all_relax_data, \(x) isTRUE(is.na(x))))))
bad_hogs <- names(which(map_lgl(all_busted_data, \(x) isTRUE(is.na(x)))))
bad_hogs

#' The errors occur as a form of this:
#'
#' ### Obtaining branch lengths and nucleotide substitution biases under the
#' nucleotide GTR model
#' Error:
#' The number of tree tips in 'XXfQsCFh.tree_0'(13) is not equal to the number
#' of sequences in the data filter associated with the tree (3).
#'
#' It's not clear to me why this is happening, so these HOGs are filtered out.
#'


# Remove these failed ones.
busted_data <- all_busted_data[map_lgl(all_busted_data, \(x) !isTRUE(is.na(x)))]
relax_data <- all_relax_data[map_lgl(all_relax_data, \(x) !isTRUE(is.na(x)))]

# These objects take up a lot of memory, so remove them from environment:
rm(all_busted_data, all_relax_data); invisible(gc())



hyphy_df <- tibble(hog = names(busted_data),
                   busted_pvals = map_dbl(busted_data, \(x) x[["test results"]][["p-value"]]),
                   relax_pvals = map_dbl(relax_data, \(x) x[["test results"]][["p-value"]]),
                   #'
                   #' From docs: "A significant result of k>1 indicates that
                   #' selection strength has been intensified along the test
                   #' branches, and a significant result of k<1 indicates that
                   #' selection strength has been relaxed" along the test
                   #' branches.
                   #'
                   relax_k = map_dbl(relax_data, \(x) {
                       n <- "relaxation or intensification parameter"
                       return(x[["test results"]][[n]])
                   }),
                   busted_disc = BY_correct(busted_pvals, hog),
                   relax_disc = BY_correct(relax_pvals, hog))



hog_focal_go_df <- "_data/hyphy-focal-hog-go.csv" |>
    read_csv(col_types = cols()) |>
    select(-offspring) |>
    rename(hog = hogs) |>
    mutate(hog = str_split(hog, ";")) |>
    unnest(hog) |>
    filter(hog %in% hyphy_df$hog) |>
    mutate(busted_disc = map_lgl(hog, \(h) hyphy_df$busted_disc[hyphy_df$hog == h]),
           relax_disc = map_lgl(hog, \(h) hyphy_df$relax_disc[hyphy_df$hog == h]))


#' Names of HOGs where null is rejected for both models.
sign_hogs <- hog_focal_go_df |>
    filter(busted_disc & relax_disc) |>
    getElement("hog") |>
    unique()


hog_focal_go_df |>
    group_by(go, term) |>
    rename(GO = go, Description = term) |>
    summarize(n = n(),
              BUSTED = mean(busted_disc) * 100,
              RELAX = mean(relax_disc) * 100,
              Both = mean(busted_disc & relax_disc) * 100,
              .groups = "drop") |>
    gt() |>
    cols_label(n = "{{*n*}}") |>
    fmt_number(BUSTED:Both, decimals = 1) |>
    tab_spanner(label = "% significant", columns = BUSTED:Both) |>
    fmt_integer(n) |>
    opt_table_font(font = "Helvetica") |>
    tab_style(style = list(cell_text(align = "center")),
              locations = cells_column_labels(-Description)) |>
    tab_options(data_row.padding = 1.5,
                column_labels.padding = 1.5,
                table.font.size = 10,
                column_labels.border.top.width = 2,
                column_labels.border.top.color = "black",
                column_labels.border.bottom.width = 1,
                column_labels.border.bottom.color = "black",
                table_body.border.top.width = 1,
                table_body.border.top.color = "black",
                table_body.border.bottom.width = 2,
                table_body.border.bottom.color = "black",
                table_body.hlines.width = 0,
                table.font.color = "black") |>
    gtsave("_figures/hyphy-go-sign-percents.pdf")

plot_crop("_figures/hyphy-go-sign-percents.pdf", quiet = TRUE)





hyphy_df |>
    filter(hog %in% sign_hogs) |>
    select(1:4) |>
    mutate(desc = case_when(hog == "N0.HOG0006267" ~ "exportin",
                            hog == "N0.HOG0001010" ~ "peroxiredoxin-1",
                            hog == "N0.HOG0005696" ~ "sulfonylurea receptor",
                            hog == "N0.HOG0007992" ~ "F-box and wd40 domain protein 7",
                            hog == "N0.HOG0007399" ~ "Serrate protein",
                            hog == "N0.HOG0007118" ~ "dynein heavy chain",
                            hog == "N0.HOG0007838" ~ "prohibitin",
                            hog == "N0.HOG0008783" ~ "protein kinase C",
                            hog == "N0.HOG0001074" ~ "ribosomal protein S6 kinase",
                            hog == "N0.HOG0004934" ~ "4-aminobutyrate aminotransferase",
                            hog == "N0.HOG0006858" ~ "Vav guanine nucleotide exchange factor",
                            hog == "N0.HOG0006910" ~ "ETS domain-containing protein",
                            hog == "N0.HOG0007233" ~ "Heat shock protein 83",
                            hog == "N0.HOG0007347" ~ "Cyclic-nucleotide-gated cation channel",
                            TRUE ~ NA)) |>
    select(hog, desc, everything()) |>
    gt() |>
    cols_label(hog = "HOG",
               desc = "Description",
               busted_pvals = "{{*P*-value}}",
               relax_pvals = "{{*P*-value}}",
               relax_k = "{{*k*}}") |>
    fmt_scientific(contains("pvals"), decimals = 2) |>
    fmt_number(relax_k, decimals = 2) |>
    tab_spanner(label = "BUSTED", columns = contains("busted")) |>
    tab_spanner(label = "RELAX", columns = contains("relax")) |>
    opt_table_font(font = "Helvetica") |>
    tab_style(style = list(cell_text(align = "center")),
              locations = cells_column_labels(-desc)) |>
    tab_options(data_row.padding = 1.5,
                column_labels.padding = 1.5,
                table.font.size = 10,
                column_labels.border.top.width = 2,
                column_labels.border.top.color = "black",
                column_labels.border.bottom.width = 1,
                column_labels.border.bottom.color = "black",
                table_body.border.top.width = 1,
                table_body.border.top.color = "black",
                table_body.border.bottom.width = 2,
                table_body.border.bottom.color = "black",
                table_body.hlines.width = 0,
                table.font.color = "black") |>
    gtsave("_figures/hyphy-sign-hogs.pdf")

plot_crop("_figures/hyphy-sign-hogs.pdf", quiet = TRUE)



