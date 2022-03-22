
library(tidyverse)
library(parallel)

options(mc.cores = max(parallel::detectCores()-2L, 1L))

read_fasta <- function(fn) {

    fasta <- read_lines(fn)

    splits <- tibble(start = which(str_detect(fasta, "^>")),
                     end = c(tail(start, -1) - 1, length(fasta))) %>%
        mutate(start = start + 1)

    seqs <- pmap_chr(splits,
                     function(start, end) str_c(fasta[start:end],
                                                collapse = ""))

    nms <- map_chr(splits$start - 1, ~ fasta[.x]) %>%
        str_remove("^>")

    names(seqs) <- nms
    class(seqs) <- "fasta"

    return(seqs)
}

print.fasta <- function(x, ...) {
    cat(sprintf("fasta file with %i sequences", length(x)))
}




ref <- read_fasta("~/_data/tany_scaffolds.fasta.gz")

# To add to each position along each scaffold to allow for plotting together:
scaff_adders <- c(0L, unname(head(nchar(ref), -1))) %>% cumsum()

# Regions of Ns
N_locs <- str_locate_all(ref, "N") %>%
    set_names(names(ref)) %>%
    map(~ .x[,1]) %>%
    imap_dfr(function(.x, .y) {
        if (length(.x) == 0) return(NULL)
        d = c(0, diff(.x))
        reg_starts <- .x[c(1L, which(d > 1))]
        reg_ends <- .x[c(which(d > 1) - 1, length(.x))]
        tibble(scaff = .y, start = reg_starts, end = reg_ends)
    }) %>%
    mutate(scaff = factor(scaff, levels = names(ref)),
           xx = scaff_adders[as.integer(scaff)],
           tog_start = start + xx,
           tog_end = end + xx,
           nbp = end - start + 1) %>%
    select(-xx)





# ============================================================*
# ============================================================*

# looking at mpileup files (1bp resolution)

# ============================================================*
# ============================================================*



samp = "MyKS-19-B_S18_bwa"

#' Add positions not present in mpileup data frame,
#' check whether each position is in a region with Ns,
#' see how close they are to N region,
#' and how close they are to end of scaffold.
process_one_mp_scaff <- function(.d) {
    scaff_idx <- which(names(ref) == .d$scaff[[1]])
    seq_size <- nchar(ref[[scaff_idx]])
    # Positions on this chromosome not present in `.d`
    zeros <- (1:seq_size)[!(1:seq_size) %in% .d$pos]
    .d <- .d %>% mutate(in_mp = TRUE)
    if (length(zeros) > 0) {
        .new_d <- tibble(scaff = .d$scaff[[1]], pos = zeros, cov = 0L,
                         tog_pos = pos + scaff_adders[scaff_idx],
                         in_mp = FALSE)
        .d <- bind_rows(.d, .new_d) %>%
            arrange(pos)
    }
    .d$end_prox <- seq_size - .d$pos
    .d$inN <- FALSE
    .d$N_prox <- NA_integer_

    nl <- N_locs %>% filter(scaff == .d$scaff[[1]])
    if (nrow(nl) == 0) return(.d)

    for (i in 1:nrow(nl)) {
        .d$inN[(nl$start[i]):(nl$end[i])] <- TRUE
    }

    .d$N_prox <- 0L
    if (nl$start[1] > 1) {
        inds <- 1:(nl$start[1] - 1L)
        .d$N_prox[inds] <- nl$start[1] - inds
    }
    if (nl$end[nrow(nl)] < seq_size) {
        inds <- (nl$end[nrow(nl)] + 1L):seq_size
        .d$N_prox[inds] <- inds - nl$end[nrow(nl)]
    }
    if (nrow(nl) > 1) {
        for (i in 2:nrow(nl)) {
            end0 <- nl$end[i-1]
            start1 <- nl$start[i]
            inds <- (end0+1):((end0 + start1) %/% 2)
            .d$N_prox[inds] <- inds - end0
            inds <- ((end0 + start1) %/% 2 + 1L):(start1-1L)
            .d$N_prox[inds] <- start1 - inds
        }
    }

    return(.d)
}


# mpileup data frame
mp_df <- read_tsv(sprintf("~/_data/%s_mpileup_3col.txt.gz", samp),
                    col_types = "cii", col_names = c("scaff", "pos", "cov")) %>%
    mutate(scaff = factor(scaff, levels = names(ref)),
           tog_pos = pos + scaff_adders[as.integer(scaff)])

# Takes ~17 sec
mp_df_list <- mp_df %>%
    split(.$scaff) %>%
    mclapply(process_one_mp_scaff)

mp_df2 <- do.call(bind_rows, mp_df_list)



splits <- list(1:65, 66:75, 76:80, 81:85)



png("~/Desktop/big_cov_plot.png", width = 20, height = 16, bg = "white", units = "in", res = 300)
grid::grid.newpage()
for (i in seq_along(splits)) {
    dd <- mp_df_list[splits[[i]]] %>%
        do.call(what = bind_rows) %>%
        mutate(col_g = ifelse(inN, "N", paste(scaff)) %>%
                   factor(levels = c("N", names(ref))))
    p <- dd %>%
        ggplot(aes(tog_pos, cov, color = col_g)) +
        geom_point(alpha = 0.1, size = 0.1) +
        theme_minimal() +
        scale_color_manual(values = c("red", rep(c("gray70", "gray40"),
                                        ceiling(length(splits[[i]]) / 2))),
                           guide = "none")
    plot(p, vp = grid::viewport(x = 0, y = 1 - (i - 1) / 4,
                                 width = 1, height = 1/4,
                                 just = c("left", "top")))
    rm(p, dd); invisible(gc())
}
dev.off()



mean(mp_df2$in_mp & mp_df2$inN)
mean(mp_df2$in_mp & !mp_df2$inN)
mean(!mp_df2$in_mp & mp_df2$inN)
mean(!mp_df2$in_mp & !mp_df2$inN)

hist(mp_df2$cov[!mp_df2$inN & mp_df2$cov < 1000])

interesting_df <- mp_df2 %>%
    filter(cov == 0 & !inN)



interesting_df %>%
    ggplot(aes(tog_pos, as.integer(in_mp), color = scaff)) +
    geom_segment(data = N_locs %>% filter(scaff %in% interesting_df$scaff),
                 aes(x = tog_start, xend = tog_end, y = -1, yend = -1),
                 color = "red", size = 1) +
    geom_point(size = 0.2) +
    scale_color_manual(values = rep(c("gray70", "gray40"),
                                    ceiling(length(ref) / 2)),
                       guide = "none") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())



# (sum(nchar(ref)) - sum(mp_df$cov > 0)) - sum(N_locs$nbp)
# hist(mp_df$cov[mp_df$cov > 0 & mp_df$cov < 600])
# hist(mp_df$cov[mp_df$cov > 600])






# ============================================================*
# ============================================================*

# looking at sliding windows

# ============================================================*
# ============================================================*



seq_type = "ill2"

seq_type <- match.arg(seq_type, c("ill", "ill2", "ont"))

if (seq_type == "ill") win_file <- "MyKS-19-B_S18_bwa_mpileup_win250.txt.gz"
if (seq_type == "ill2") win_file <- "Lys-19_S14_bwa_mpileup_win250.txt.gz"
if (seq_type == "ont") win_file <- "ont_mpileup_win250.txt.gz"


win_df <- read_tsv(paste0("~/_data/", win_file), col_types = "ciiiidd") %>%
    mutate(mid = as.integer((start + end) / 2)) %>%
    mutate(scaff = factor(scaff, levels = names(ref)),
           tog_mid = mid + scaff_adders[as.integer(scaff)])



for (zoom in c(TRUE, FALSE)) {
    p <- win_df %>%
        ggplot(aes(tog_mid, median, color = scaff)) +
        geom_segment(data = N_locs,
                     aes(x = tog_start, xend = tog_end, y = -10, yend = -10),
                     color = "red", size = 1) +
        geom_point(size = 0.2) +
        scale_color_manual(values = rep(c("gray70", "gray40"),
                                        ceiling(length(ref) / 2)),
                           guide = "none") +
        theme_minimal() +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    if (zoom) p <- p + coord_cartesian(ylim = c(-20, 500))
    fn <- sprintf("~/Desktop/cov_%s%s.png", seq_type, ifelse(zoom, "_zoom", ""))
    ggsave(fn, p, width = 10, height = 4.15, bg = "white")
    cat(sprintf("Saved %s\n", fn))
}







# Number of N regions increases linearly with scaffold size:

imap_dfr(ref,
         function(.x, .y) {
             tibble(scaff = .y, nregs = sum(N_locs$scaff == .y),
                    size = nchar(ref[[.y]]))
         }) %>%
    ggplot(aes(size / 1e6, nregs)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    xlab("Size of scaffold (Mb)") +
    ylab("Number of regions with 'N'") +
    theme_minimal()


