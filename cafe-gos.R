


library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
# ^^ make sure these are loaded before tidyverse ^^

library(rrvgo)
library(org.Dm.eg.db)

source("_scripts/00-preamble.R")


#' #' Show node names of cafe tree:
#' ```r
#' library(ape)
#' library(ggtree)
#' cafe_tr <- paste0(cafe_dir, "/cafe_k8_run1/Gamma_asr.tre") |>
#'     read_lines() |>
#'     keep(~ startsWith(.x, "  TREE")) |>
#'     getElement(1) |>
#'     str_remove(".*= ") |>
#'     (\(x) read.tree(text = x))()
#' ggtree(cafe_tr) +
#'     geom_tiplab() +
#'     geom_nodelab()
#' ```
#' #' from this, I can see that Chironomidae starts at node 17 (from node 20)

.node <- rlang::sym("<17>")


#' Gene count changes:
gc_df <- paste0(cafe_dir, "/cafe_k8_run1/Gamma_change.tab") |>
    read_tsv(col_types = paste0(c("c", rep("d", 25)), collapse = "")) |>
    # change in count from node 20 to 17 (i.e., to chironomidae):
    mutate(chir_d = !!.node, hog = FamilyID) |>
    select(hog, chir_d)


#' Start by reading branch probabilities, to get p-value (at node 17),
#' then combine with `gc_df` to add change in gene count (from node 17 to 20)
#' for each HOG.
#' I had to do the column names this way because it has a weird extra
#' tab at the end of the line.
bp_file <- paste0(cafe_dir, "/cafe_k8_run1/Gamma_branch_probabilities.tab")
hog_pd_df <- read_tsv(bp_file,
                         na = "N/A", col_types =
                             paste0(c("c", rep("d", 25),"c"), collapse = ""),
                      skip = 1,
                      col_names = str_split(read_lines(bp_file)[1], "\t")[[1]] |>
                          head(-1)) |>
    rename(hog = `#FamilyID`, pval = !!.node) |>
    select(hog, pval) |>
    left_join(gc_df, by = "hog")

#' Relationship between change in gene count and p-value:
hog_pd_df |>
    ggplot(aes(abs(chir_d), -log(pval))) +
    geom_hline(yintercept = -log(0.05)) +
    geom_jitter(height = 0, width = 0.05, alpha = 0.5, shape = 1)



# GO terms for all HOGs:
hog_gos <- paste0(orthofinder_extr_dir, "/All_HOG_GO/N0-GO-by-HOG.tsv") |>
    read_tsv(col_types = cols())



#' All N0 HOGs with GO terms and descriptions.
#' No filtering for whether they rapidly evolved at node 17 yet!
chir_hogs <- hog_pd_df |>
    mutate(go = map_chr(hog, \(h) hog_gos$go[hog_gos$hog == h]) |> toupper()) |>
    filter(!is.na(go)) |>
    mutate(go = str_split(go, ";")) |>
    unnest(go) |>
    mutate(term = suppressMessages(AnnotationDbi::select(GO.db, go, "TERM")[,"TERM"])) |>
    # remove GO terms that weren't found in DB bc they're obsolete:
    filter(!is.na(term)) |>
    mutate(phred = -log(pval))





# =============================================================================*
# =============================================================================*
# Over-representation test:
# =============================================================================*
# =============================================================================*



term2gene <- chir_hogs |>
    select(go, hog)
term2name <- chir_hogs |>
    select(go, term)

genes <- chir_hogs |>
    distinct(hog, pval, chir_d) |>
    arrange(pval) |>
    filter(pval < 0.001, chir_d > 0) |>
    getElement("hog") |>
    # base::`[`(1:20) |>
    identity()




overrep <- enricher(genes,
                    pvalueCutoff = 0.1, pAdjustMethod = "BH",
                    minGSSize = 1,
                    TERM2GENE=term2gene, TERM2NAME=term2name)
overrep |> as_tibble()


or_df <- overrep |>
    as_tibble() |>
    #' Ontology:
    #'   * MF - molecular function
    #'   * BP - biological process
    #'   * CC - cellular component
    mutate(ont = Ontology(ID))

or_bp_scores <- or_df |>
    filter(ont == "BP") |>
    mutate(geneID = geneID |> str_split("/")) |>
    unnest(geneID) |>
    mutate(d_chr = map_int(geneID, \(x) hog_pd_df$chir_d[hog_pd_df$hog == x])) |>
    group_by(ID) |>
    summarize(d_chr = sum(d_chr)) |>
    (\(x) {z <- (x$d_chr); names(z) <- x$ID; return(z)})()

#' `org.Dm.eg.db` is the genome wide annotation for *Drosophila melanogaster*
or_bp_sim_mat <- calculateSimMatrix(names(or_bp_scores),
                                    "org.Dm.eg.db", ont = "BP")

or_bp_red <- reduceSimMatrix(or_bp_sim_mat,
                          scores = or_bp_scores,
                          orgdb = "org.Dm.eg.db") |>
    as_tibble()
or_bp_red

treemapPlot(or_bp_red)
# scatterPlot(or_bp_sim_mat, or_bp_red, onlyParents = TRUE)
wordcloudPlot(or_bp_red, min.freq=1, colors="black")





