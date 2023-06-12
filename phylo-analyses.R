library(tidyverse)
library(ape)
library(phyr)
library(treeio)
library(patchwork)

theme_set(theme_classic())

save_plot <- function(n, p, w, h, ...) {
    cairo_pdf(sprintf("~/Stanford_Drive/UW/midgenomes/%s.pdf", n),
              width = w, height = h, bg = NA, ...)
    plot(p)
    dev.off()
}


#' Chironomid assembly data:
chir_df <-
"Species\tTotal length (Mb)\tProtein-coding genes\tRepeat content (%)\tMean intron length (bp)
Tanytarsus gracilentus\t91.83\t15700\t8.01\t451.6
Chironomus tepperi\t236.64\tNA\tNA\tNA
Chironomus riparius\t178.17\t13449\t9.14\t508.4
Chironomus tentans\t213.46\t15120\t10\t1103
Polypedilum vanderplanki\t118.97\t18990\t5.01\t532.6
Polypedilum pembai\t122.92\t15068\tNA\tNA
Belgica antarctica\t89.58\t13517\t0.4\t333
Clunio marinus\t85.49\t22620\t6.86\t383.5
Propsilocerus akamusi\t85.84\t11942\t6.38\t528.6
Parochlus steinenii\t143.57\t14838\t12.38\t567.1" |>
    read_tsv(col_types = "cdddd") |>
    set_names(c("species", "length", "genes", "repeats", "intron_len")) |>
    filter(!is.na(genes) | !is.na(repeats) | !is.na(intron_len))

#' Version without Polypedilum pembai bc it doesn't have repeats or intron data
chir_df2 <- chir_df |>
    filter(!is.na(repeats) & !is.na(intron_len))

#' A third one also without Chironomus tentans to see how much it affects
#' the intron analysis bc it has such a high value.
chir_df3 <- chir_df2 |>
    filter(species != "Chironomus tentans")


#' time tree:
chir_tr <- read.mcmctree("~/_data/phylo/chir_mcmctree/mcmc_1/FigTree.tre") |>
    .@phylo
chir_tr$tip.label <- gsub("_", " ", chir_tr$tip.label)
chir_tr <- keep.tip(chir_tr, chir_df$species)

chir_tr2 <- keep.tip(chir_tr, chir_df2$species)
chir_tr3 <- keep.tip(chir_tr, chir_df3$species)



m_lg <- cor_phylo(~ length + genes, ~ species, chir_tr, data = chir_df,
                  constrain_d = TRUE)
m_lg0 <- update(m_lg, no_corr = TRUE)
m_lr <- cor_phylo(~ length + repeats, ~ species, chir_tr2, data = chir_df2,
                  constrain_d = TRUE)
m_lr0 <- update(m_lr, no_corr = TRUE)
m_li <- cor_phylo(~ length + intron_len, ~ species, chir_tr2, data = chir_df2,
                  constrain_d = TRUE)
m_li0 <- update(m_li, no_corr = TRUE)
# Trying last one without Chironomus tentans:
m_li2 <- cor_phylo(~ length + intron_len, ~ species, chir_tr3,
                  constrain_d = TRUE,
                  data = chir_df3)
m_li20 <- update(m_li2, no_corr = TRUE)

(lrt_lg <- -2 * (m_lg0$logLik - m_lg$logLik))
pchisq(lrt_lg, df = 1, lower.tail = FALSE)

(lrt_lr <- -2 * (m_lr0$logLik - m_lr$logLik))
pchisq(lrt_lr, df = 1, lower.tail = FALSE)

(lrt_li <- -2 * (m_li0$logLik - m_li$logLik))
pchisq(lrt_li, df = 1, lower.tail = FALSE)

(lrt_li2 <- -2 * (m_li20$logLik - m_li2$logLik))
pchisq(lrt_li2, df = 1, lower.tail = FALSE)


m_lg$d
m_lr$d
m_li$d



lg_p <- chir_df |>
    ggplot(aes(length, genes / 1000)) +
    geom_point() +
    xlab(NULL) +
    ylab("Thousands of genes")

lr_p <- chir_df |>
    ggplot(aes(length, repeats)) +
    geom_point(na.rm = TRUE) +
    xlab("Genome size (Mb)") +
    ylab("Repeats (%)")

li_p <- chir_df |>
    ggplot(aes(length, intron_len)) +
    geom_point(na.rm = TRUE) +
    xlab(NULL) +
    ylab("Mean intron length (bp)")


length_p <- lg_p +
    lr_p +
    li_p +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"))

save_plot(n = "size_structure", p = length_p, w = 6.5, h = 2.5)





library(phylolm)

pl_chir_df <- chir_df |> as.data.frame()
rownames(pl_chir_df) <- pl_chir_df$species

pl_chir_df2 <- chir_df2 |> as.data.frame()
rownames(pl_chir_df2) <- pl_chir_df2$species

phylolm(genes ~ length, pl_chir_df, chir_tr, model = "OUfixedRoot") |> summary()
phylolm(repeats ~ length, pl_chir_df2, chir_tr2) |> summary()
phylolm(intron_len ~ length, pl_chir_df2, chir_tr2) |> summary()



