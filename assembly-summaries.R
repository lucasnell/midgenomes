
#'
#' Create figures to summarize assembly stats.
#'

library(tidyverse)

theme_set(theme_classic() +
              theme(strip.background = element_blank(),
                    axis.text = element_text(color = "black")))


save_plot <- function(n, p, w, h, ...) {
    cairo_pdf(sprintf("~/Box Sync/midgenomics/genome-report/%s.pdf", n),
              width = w, height = h, bg = NA, ...)
    plot(p)
    dev.off()
}

# Adjust data frame to use for plotting BUSCO info.
prep_busco <- function(.df) {
    .df %>%
        mutate(across(all_of(c("BUSCO_C-S", "BUSCO_C-D", "BUSCO_F", "BUSCO_M")),
                      ~ 100 * .x / BUSCO_n)) %>%
        select(-busco_text, -BUSCO_C, -BUSCO_n) %>%
        pivot_longer(starts_with("BUSCO"), names_to = "b_class") %>%
        mutate(b_class = case_when(b_class == "BUSCO_C-S" ~ "complete (C) and single-copy (S)",
                                   b_class == "BUSCO_C-D" ~ "complete (C) and duplicated (D)",
                                   b_class == "BUSCO_F" ~ "fragmented (F)",
                                   b_class == "BUSCO_M" ~ "missing (M)",
                                   TRUE ~ NA_character_) %>%
                   factor(levels = c("complete (C) and single-copy (S)",
                                     "complete (C) and duplicated (D)",
                                     "fragmented (F)", "missing (M)") %>%
                              rev()))
}
# Fill scale for BUSCO classes
busco_fill_scale <- scale_fill_manual(
    NULL,
    values = c("dodgerblue", "dodgerblue4", "gold", "firebrick"),
    limits = c("complete (C) and single-copy (S)",
               "complete (C) and duplicated (D)",
               "fragmented (F)", "missing (M)"))


# ============================================================================*
# ============================================================================*

# main assemblies ----

# ============================================================================*
# ============================================================================*



#' Below, I've directly inserted comma-delimited text representing a table
#' of summary stats on the various assemblies.
#'
#' Column descriptions:
#' - assembler: Program used for the initial assembly, where `quickmerge`
#'   indicates the best permutation from the `quickmerge` summary below.
#' - step: Program used in each step, where ">" indicates that this
#'   program was run on the output from the previous step
#' - size: Total size of assembly (Mb)
#' - n_seqs: Number of contigs
#' - N50: Contig N50 (Mb)
#' - min: Minimum contig size (bp)
#' - max: Max contig size (Mb)
#' - total_N: Total N present in assembly (bp)
#'   Because this is 0 for all assemblies, it's not show in figures.
#' - BUSCO_C: Number of complete BUSCO genes
#' - BUSCO_C-S: Number of complete & single-copy BUSCO genes
#' - BUSCO_C-D: Number of complete & duplicated BUSCO genes
#' - BUSCO_F: Number of fragmented BUSCO genes
#' - BUSCO_M: Number of missing BUSCO genes
#' - BUSCO_n: Total number of BUSCO genes
#'

summ_df <-
"assembler,step,size,n_seqs,N50,min,max,total_N,BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n
NECAT,NECAT,93.488293,198,2.154540,619,3.883527,0,2912,2845,67,88,285,3285
NECAT,> Racon,93.598262,192,2.161385,1062,3.895083,0,2933,2867,66,86,266,3285
NECAT,> medaka,93.848166,192,2.165007,1070,3.907711,0,2963,2894,69,60,262,3285
NECAT,> purge_dups,91.895536,121,2.125084,2004,3.907711,0,2963,2913,50,59,263,3285
SMARTdenovo,SMARTdenovo,95.110196,159,1.396411,7167,5.335390,0,2880,2798,82,116,289,3285
SMARTdenovo,> Racon,95.422357,159,1.406639,5040,5.373509,0,2949,2866,83,87,249,3285
SMARTdenovo,> medaka,95.468663,159,1.411120,5077,5.366456,0,2988,2898,90,58,239,3285
SMARTdenovo,> purge_dups,92.841149,119,1.519394,6612,5.366456,0,2985,2930,55,59,241,3285
NextDenovo,NextDenovo,93.033974,88,2.230847,16933,6.769189,0,2945,2906,39,77,263,3285
NextDenovo,> Racon,92.804737,88,2.230921,16037,6.761443,0,2950,2917,33,81,254,3285
NextDenovo,> medaka,92.612204,88,2.228801,16049,6.760912,0,2987,2955,32,57,241,3285
Flye,Flye,97.111699,1114,0.927999,495,3.695663,0,2946,2860,86,80,259,3285
Flye,> Racon,97.507934,1052,0.936236,92,3.724113,0,2947,2868,79,80,258,3285
Flye,> medaka,97.873717,1052,0.939342,92,3.738528,0,2987,2907,80,59,239,3285
Flye,> purge_dups,91.893637,333,0.993865,203,3.738528,0,2987,2946,41,59,239,3285
quickmerge,quickmerge,92.331558,45,7.052781,16049,12.427199,0,2991,2959,32,57,237,3285
quickmerge,> medaka,92.029278,45,7.024462,16040,12.399045,0,2993,2963,30,57,235,3285
quickmerge,> NextPolish,91.827299,45,7.014169,16039,12.374712,0,3009,2974,35,49,227,3285" %>%
    read_csv(col_types = "ccdididiiiiiii") %>%
    # mutate(step = ifelse(step %in% assembler, "assemble", step) %>%
    #            str_remove("^> ") %>%
    #            # factor(levels = c("assemble", "Racon", "medaka", "purge_dups")),
    #            factor(levels = c("purge_dups", "medaka", "Racon", "assemble")),
    #        busco_text = sprintf("C:%i [S:%i, D:%i], F:%i, M:%i, n:%i",
    #                             BUSCO_C,`BUSCO_C-S`,`BUSCO_C-D`,BUSCO_F,
    #                             BUSCO_M,BUSCO_n))
    mutate(busco_text = sprintf("C:%i [S:%i, D:%i], F:%i, M:%i, n:%i",
                                BUSCO_C,`BUSCO_C-S`,`BUSCO_C-D`,BUSCO_F,
                                BUSCO_M,BUSCO_n),
           step_id = ifelse(step == assembler, step, paste(assembler, step)),
           step_id = factor(step_id, levels = rev(step_id)),
           assembler = factor(assembler, levels = unique(assembler)))


busco_labeller <- function(x, .width = 15) {
    repl <- gsub("NECAT|SMARTdenovo|NextDenovo|Flye|quickmerge", " ", x)
    newx <- ifelse(!str_detect(x, "\\>"), x, repl)
    str_pad(newx, .width, side = "right")
}

busco_p <- summ_df %>%
    prep_busco() %>%
    ggplot(aes(step_id, value)) +
    geom_col(aes(fill = b_class)) +
    geom_text(data = summ_df,
              aes(y = 5, label = busco_text),
              hjust = 0, size = 6 / 2.81) +
    busco_fill_scale +
    scale_x_discrete(labels = busco_labeller) +
    ylab("% BUSCOs") +
    coord_flip() +
    facet_wrap(~ assembler, ncol = 1, scales = "free_y") +
    guides(fill = guide_legend(nrow = 2)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(family = "Courier", hjust = 0, size = 8),
          legend.position = "none",
          panel.spacing.y = unit(0, "lines"),
          strip.text = element_blank()) +
    NULL

save_plot("busco", busco_p, 5, 4.5)


summ_df %>%
    select(step, size, n_seqs, N50, min, max) %>%
    mutate(across(where(is.double), ~ sprintf("%.2f", .x))) %>%
    format_csv() %>%
    cat()




# ============================================================================*
# ============================================================================*

# quickmerge permutations ----

# ============================================================================*
# ============================================================================*


#'
#' Same as above, except for the outputs from `quickmerge`.
#'
#' The only structural difference is that instead of columns `assembler` and
#' `step`, there is the `perm` column that indicates the permutation
#' (because order matters in `quickmerge`) of assemblies used.
#'

merge_summ_df <-
"perm,size,n_seqs,N50,min,max,total_N,BUSCO_C,BUSCO_C-S,BUSCO_C-D,BUSCO_F,BUSCO_M,BUSCO_n
flye + necat,91.446177,180,6.049977,203,13.960534,0,2955,2912,43,59,271,3285
flye + next,93.316631,156,6.498738,203,18.171453,0,2970,2928,42,55,260,3285
flye + smart,89.868646,192,4.503389,203,11.068235,0,2897,2851,46,53,335,3285
necat + flye,91.451904,103,2.757536,2004,13.921865,0,2960,2912,48,57,268,3285
necat + next,93.765062,106,2.772439,2004,6.760912,0,2958,2848,110,58,269,3285
necat + smart,92.408460,96,3.069729,2004,13.296413,0,2966,2899,67,56,263,3285
next + flye,92.344005,68,3.305925,16049,12.427199,0,2990,2958,32,57,238,3285
next + necat,93.491206,75,2.889028,17544,6.760912,0,2959,2879,80,58,268,3285
next + smart,92.448745,62,3.305925,17544,9.013862,0,2988,2958,30,58,239,3285
smart + flye,92.332293,69,6.302447,6612,11.252327,0,2985,2928,57,58,242,3285
smart + necat,92.473987,53,5.307899,6612,13.239600,0,2941,2858,83,57,287,3285
smart + next,92.111530,56,5.384500,6612,11.627627,0,2987,2954,33,56,242,3285
(next + flye) + (next + smart),94.800438,60,6.052544,16049,12.427199,0,2990,2878,112,57,238,3285
(next + smart) + (next + flye),92.287285,54,6.052544,17544,12.427199,0,2990,2959,31,57,238,3285
(next + flye) + (smart + necat),93.060578,40,7.052781,17544,13.245146,0,2977,2902,75,59,249,3285
(smart + necat) + (next + flye),90.208033,41,8.576887,6612,14.248256,0,2907,2854,53,56,322,3285
(next + flye) + (smart + next),92.331558,45,7.052781,16049,12.427199,0,2991,2959,32,57,237,3285
(smart + next) + (next + flye),92.008947,46,7.052781,6612,12.428575,0,2988,2953,35,55,242,3285" %>%
    read_csv(col_types = "cdididiiiiiii") %>%
    mutate(busco_text = sprintf("C:%i [S:%i, D:%i], F:%i, M:%i, n:%i",
                                BUSCO_C,`BUSCO_C-S`,`BUSCO_C-D`,BUSCO_F,
                                BUSCO_M,BUSCO_n))

#' Versions used downstream for each round.
#' All but the last item are the permutations used for second round.
#' The last item shows the second-round permutation used to produce
#' final assembly.
besties <- c("next + flye", "next + smart", "smart + necat", "smart + next",
             "(next + flye) + (smart + next)")
merge_summ_df <- merge_summ_df %>%
    mutate(perm = ifelse(perm %in% besties, paste("*", perm), paste(" ", perm)),
           # Change all to 4 characters:
           perm = str_replace_all(perm, "necat", "neca") %>%
               str_replace_all("smart", "smar"),
           # Order in CSV is reverse from what should be plotted:
           perm = factor(perm, levels = rev(perm)))


merge_busco_p <- merge_summ_df %>%
    prep_busco() %>%
    ggplot(aes(perm, value)) +
    geom_col(aes(fill = b_class)) +
    geom_text(data = merge_summ_df, aes(y = 5, label = busco_text), hjust = 0,
              size = 6 / 2.81) +
    busco_fill_scale +
    ylab("% BUSCOs") +
    coord_flip() +
    guides(fill = guide_legend(nrow = 2)) +
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(family = "Courier")) +
    NULL





