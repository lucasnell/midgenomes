
Rcpp::sourceCpp("helpers.cpp")


# From decimal degrees to UTM, assuming it's Iceland and using WGS84
to_utm <- function(.df, .lat = "lat", .lon = "lon") {
    .cord.dec <- sp::SpatialPoints(cbind(.df[[.lon]], .df[[.lat]]),
                                   proj4string = sp::CRS("+proj=longlat"))
    .cord.UTM <- sp::spTransform(.cord.dec,
                                 sp::CRS("+proj=utm +zone=28 ellps=WGS84"))
    .df[[.lon]] <- .cord.UTM@coords[,1]
    .df[[.lat]] <- .cord.UTM@coords[,2]
    return(.df)
}

#' This function does two things to the input `freq_df` data frame:
#'
#' (1) Filter for biallelic loci only.
#' (2) Correct allele frequencies to make sure they refer to the same allele
#'     for all pools.
#'
correct_biallelic_freqs <- function(freq_df, count_df) {

    .pools <- colnames(freq_df)[-1:-3]

    # Allele counts
    if (is.character(count_df)) {
        # Takes ~9 sec
        count_df <- read_tsv(count_df,
                             col_names = c("scaff", "pos", "ref", .pools),
                             col_types = paste(c("cic", rep("c", length(.pools))),
                                               collapse = ""),
                             progress = FALSE) %>%
            mutate(across(all_of(.pools), split_sync_strings)) %>%
            identity()
    } else if (!inherits(count_df, "data.frame")) {
        stop("count_df must be string or data frame")
    }

    stopifnot(all(colnames(freq_df) == colnames(count_df)))

    n_alleles <- count_alleles(count_df %>% select(all_of(.pools)) %>% as.list())

    # Counts and frequencies of biallelic alleles:
    allele_freqs <- freq_df %>%
        filter(n_alleles == 2) %>%
        select(all_of(.pools)) %>%
        as.matrix()
    allele_counts <- count_df %>%
        filter(n_alleles == 2) %>%
        select(all_of(.pools)) %>%
        as.list()

    # allele each allele frequency refers to:
    which_alleles <- which_allele(allele_counts, allele_freqs)
    # alleles whose frequencies are >0 for each locus
    nz_alleles <- nonzero_alleles(allele_counts)

    # This should always be zero:
    nz_which_mismatch <- map_lgl(1:nrow(nz_alleles),
                                 function(i) !all(which_alleles[i,] %in% nz_alleles[i,])) %>%
        sum()
    stopifnot(nz_which_mismatch == 0)

    # Because of the above check, we can simply use 1-p when p refers
    # to the wrong allele ("wrong" = refers to a different allele from what
    # most of the others refer to)

    mode_avg <- function(v) {
        uniq_v <- unique(v)
        uniq_v[which.max(tabulate(match(v, uniq_v)))]
    }

    to_switch <- map(1:nrow(nz_alleles),
                     function(i) {
                         m <- mode_avg(which_alleles[i,])
                         which_alleles[i,] != m
                     }) %>%
        do.call(what = rbind)
    cat(sprintf("Fixing %s (%s %%) frequencies\n",
                prettyNum(sum(to_switch), big.mark=","),
                format(100 * mean(to_switch), digits = 3)))
    allele_freqs2 <- allele_freqs
    allele_freqs2[to_switch] <- 1 - allele_freqs[to_switch]

    # Now convert back to tibble for output:
    out <- freq_df %>%
        filter(n_alleles == 2)
    # Fill corrected allele frequencies:
    for (p in .pools) out[[p]] <- allele_freqs2[,p]

    return(out)

}

