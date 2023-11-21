
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


#'
#' The same 12 HOGs failed for both tests:
#'
#' - N0.HOG0005474
#' - N0.HOG0005567
#' - N0.HOG0006093
#' - N0.HOG0006291
#' - N0.HOG0006401
#' - N0.HOG0006823
#' - N0.HOG0006968
#' - N0.HOG0007232
#' - N0.HOG0007353
#' - N0.HOG0007393
#' - N0.HOG0007749
#' - N0.HOG0008307
#'
#'
#' The errors occur as a form of this:
#'
#' ### Obtaining branch lengths and nucleotide substitution biases under the
#' nucleotide GTR model
#' Error:
#' The number of tree tips in 'XXfQsCFh.tree_0'(13) is not equal to the number
#' of sequences in the data filter associated with the tree (3).
#'
#' This happens because they "failed to align to any of the in-frame references"
#' in the first step of the codon-aware alignment
#' (`hyphy /opt/codon-msa/pre-msa.bf --input HOG-NAME.fasta`)
#' that happens inside `hyphy-align.sh`.
#' This results in alignments with fewer species than tips in the phylogeny.
#' These HOGs are discarded from our analysis.
#'


# Read HyPhy results and return NA for failed HOGs:
read_hyphy <- function(x, .test) {
    tryCatch(read_json(paste0(dirs[[paste0("hyphy_",tolower(.test))]], "/", x)),
             error = function(e) NA)
}

busted_data <- list.files(paste0(dirs$hyphy_busted, "/", ""), "*.json")|>
    set_names(\(x) str_remove(x, ".json$")) |>
    map(read_hyphy, .test = "busted") |>
    discard(\(x) isTRUE(is.na(x)))

relax_data <- list.files(paste0(dirs$hyphy_relax, "/", ""), "*.json") |>
    set_names(\(x) str_remove(x, ".json$")) |>
    map(read_hyphy, .test = "relax") |>
    discard(\(x) isTRUE(is.na(x)))





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

hog_focal_go_p_df <- hog_focal_go_df |>
    group_by(go, term) |>
    summarize(BUSTED = mean(busted_disc),
              RELAX = mean(relax_disc),
              both = mean(busted_disc & relax_disc),
              .groups = "drop") |>
    pivot_longer(BUSTED:RELAX) |>
    mutate(name = factor(name, levels = rev(c("BUSTED", "Both", "RELAX"))),
           term = sprintf("%s (%s)", term, go)) |>
    arrange(go, term, name) |>
    group_by(go, term) |>
    mutate(start = c(0, value[1] - both[1]),
           end = c(value[1], value[1] + value[2] - both[2])) |>
    ungroup() |>
    mutate(term = factor(term, levels = sort(unique(term))),
           term_int = as.integer(term))

ho_focal_go_lab_df <- hog_focal_go_df |>
    group_by(go, term) |>
    summarize(n = n(),
              value = mean(busted_disc | relax_disc) + 0.03,
              .groups = "drop") |>
    mutate(term = sprintf("%s (%s)", term, go),
           term = factor(term, levels = sort(unique(term))),
           term_int = as.integer(term))


hog_focal_go_p <- hog_focal_go_p_df |>
    mutate(term_int = case_when(str_detect(term, "cold|anoxia|ionizing") ~ term_int,
                                name == "RELAX" ~ term_int + 0.2,
                                TRUE ~ term_int - 0.2)) |>
    ggplot(aes(y = term_int))  +
    geom_segment(aes(yend = term_int, x = start, xend = end, color = name),
                 linewidth = 3) +
    geom_text(data = ho_focal_go_lab_df,
              aes(x = value, label = n), hjust = 0, size = 8/2.8) +
    scale_x_continuous("Significant HOGs", labels = scales::percent,
                       limits = c(0, 1.1), breaks = 0.25* 0:4) +
    scale_y_continuous(NULL, breaks = 1:length(levels(hog_focal_go_p_df$term)),
                       labels = levels(hog_focal_go_p_df$term),
                       limits = c(1, length(levels(hog_focal_go_p_df$term))) +
                           c(-1, 1) * 0.25) +
    scale_color_manual(values = viridisLite::magma(2, begin = 0.3, end = 0.8)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.position = "none",
          axis.text.y = element_text(color = "black"))

save_plot("hyphy-focal-go-percent", hog_focal_go_p, w = 4, h = 1.6, .png = FALSE)







#' Names of HOGs where null is rejected for both models.
sign_hogs <- hog_focal_go_df |>
    filter(busted_disc & relax_disc) |>
    getElement("hog") |>
    unique()

#'
#' Below, I replace zeros with <1e-26 because this is more accureate as
#' indicated by one of the authors:
#' https://github.com/veg/hyphy/issues/988#issuecomment-504519957
#'
hyphy_sign_hogs_table <- hyphy_df |>
    filter(hog %in% sign_hogs) |>
    select(1:4) |>
    mutate(desc = case_when(hog == "N0.HOG0006267" ~ "exportin",
                            hog == "N0.HOG0001010" ~ "peroxiredoxin-1",
                            hog == "N0.HOG0005696" ~ "sulfonylurea receptor",
                            hog == "N0.HOG0007992" ~ "F-box/WD-40 domain protein 7",
                            hog == "N0.HOG0007399" ~ "SERRATE protein",
                            hog == "N0.HOG0007118" ~ "dynein heavy chain",
                            hog == "N0.HOG0007838" ~ "prohibitin",
                            hog == "N0.HOG0008783" ~ "protein kinase C",
                            hog == "N0.HOG0001074" ~ "ribosomal protein S6 kinase",
                            hog == "N0.HOG0004934" ~ "4-aminobutyrate aminotransferase",
                            hog == "N0.HOG0006858" ~ "Vav guanine nucleotide exchange factor",
                            hog == "N0.HOG0006910" ~ "ETS domain-containing protein",
                            hog == "N0.HOG0007233" ~ "heat shock protein 83",
                            hog == "N0.HOG0007347" ~ "cyclic-nucleotide-gated cation channel",
                            TRUE ~ NA)) |>
    # So that sub_small_vals works below:
    mutate(across(contains("pvals"), \(x) ifelse(x == 0, 1e-14, x))) |>
    select(hog, desc, everything()) |>
    arrange(relax_k) |>
    gt() |>
    cols_label(hog = "HOG",
               desc = "Description",
               busted_pvals = "{{*P*-value}}",
               relax_pvals = "{{*P*-value}}",
               relax_k = "{{*k*}}") |>
    fmt_scientific(contains("pvals"), decimals = 2) |>
    sub_small_vals(columns = contains("pvals"),
                   threshold = 1e-13,
                   small_pattern = md("<10<sup>−26</sup>")) |>
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
                table.font.color = "black")

hyphy_sign_hogs_table |>
    gtsave("_figures/hyphy-sign-hogs.pdf")

plot_crop("_figures/hyphy-sign-hogs.pdf", quiet = TRUE)



