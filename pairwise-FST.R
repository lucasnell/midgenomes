library(tidyverse)

theme_set(theme_classic())
theme_update(axis.text = element_blank(),
             axis.title = element_blank(),
             axis.ticks = element_blank(),
             axis.line = element_line(size = 1.5))

f <- function(.x) {
    x_sin <- sin((.x / diff(range(.x)) - 0.5/3) * 3 * pi) + 1
    x <- .x / 10
    flr <- (x %/% 2) * 2
    x <- x - flr
    nc <- x * pi
    z <- sin(nc)
    z <- z + x_sin * 0.5
    return(z)
}

# curve(f(x), 0, 100)

pred_prey_df <- tibble(x = seq(0, 100, length.out = 1001),
                       y = f(x),
                       y2 = f(x - 5))
bal_sel_df <- tibble(x = seq(0, 100, length.out = 1001),
       y = f(x - 10))

pred_prey_p <- pred_prey_df %>%
    ggplot(aes(x, y)) +
    geom_line(size = 1, color = "dodgerblue3") +
    geom_line(aes(y = y2), size = 1, color = "gold")

bal_sel_p <- bal_sel_df %>%
    ggplot(aes(x, y)) +
    geom_line(size = 1, color = "firebrick")

ggsave("~/Desktop/pred_prey.pdf", pred_prey_p, width = 5, height = 2.5)
ggsave("~/Desktop/bal_sel.pdf", bal_sel_p, width = 5, height = 2.5)





library(tidyverse)
library(readxl)
library(rgdal)
library(broom)
library(poolfstat)
library(vegan)
library(rstan)

options(mc.cores = max(parallel::detectCores()-2, 1))
rstan_options(auto_write = TRUE)

# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")


# From decimal degrees to UTM, assuming it's Iceland and using WGS84
to_utm <- function(.df, .lat = "lat", .lon = "lon") {
    .cord.dec <- SpatialPoints(cbind(.df[[.lon]], .df[[.lat]]),
                               proj4string=CRS("+proj=longlat"))
    .cord.UTM <- spTransform(.cord.dec, CRS("+proj=utm +zone=28 ellps=WGS84"))
    .df[[.lon]] <- .cord.UTM@coords[,1]
    .df[[.lat]] <- .cord.UTM@coords[,2]
    return(.df)
}

samp_df <- read_csv("~/Box Sync/midgenomics/full-DNA-info.csv",
                    col_types = "cfcddidcccddidiDccldd") %>%
    # Convert to more accurate site names:
    mutate(site = case_when(site == "SN" ~ "Syðri Neslönd",
                            site == "KS" ~ "Kálfaströnd",
                            site == "H" ~ "Haganes",
                            site == "V" ~ "Vindbelgur",
                            site == "BR" ~ "Fellshóll",
                            TRUE ~ NA_character_),
           location = ifelse(lake == "Mývatn", site, lake)) %>%
    filter(read == 1) %>%
    filter(date >= as.Date("2019-01-01")) %>%
    select(biotech_id, n_adults, date, lake, site, location, lat, lon) %>%
    to_utm()


pools <- read_lines("~/_data/snape/space_sync_masked.names.gz") %>%
    map_chr(~ samp_df$location[samp_df$biotech_id == .x])
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$location == .x][1])

# Takes ~16 sec with 6 threads
pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/space_sync_masked_noblanks.sync.gz",
                             poolsizes = n_adults * 2,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = getOption("mc.cores", 2L))


# Identity method seems to work better for dealing with low number of
# individuals in KS pool
pw_fst <- compute.pairwiseFST(pool_dat, method = "Identity",
                              nsnp.per.bjack.block = 0)
# nsnp.per.bjack.block = pool_dat@nsnp %/% 500)

max_fst <- max(pw_fst@PairwiseFSTmatrix, na.rm = TRUE)
max_inds <- which(pw_fst@PairwiseFSTmatrix == max_fst)
max_cols <- ceiling(max_inds / nrow(pw_fst@PairwiseFSTmatrix))
max_rows <- max_inds %% nrow(pw_fst@PairwiseFSTmatrix)
pw_fst@PairwiseFSTmatrix[max_rows, max_cols]

# Below is missing a negative Fst, not sure why.
# neg_inds <- which(pw_fst@PairwiseFSTmatrix < 0)
# neg_cols <- unique(ceiling(neg_inds / nrow(pw_fst@PairwiseFSTmatrix)))
# neg_rows <- unique(neg_inds %% nrow(pw_fst@PairwiseFSTmatrix))
# pw_fst@PairwiseFSTmatrix[neg_rows, neg_cols]
# sum(pw_fst@PairwiseFSTmatrix[neg_rows, neg_cols] < 0, na.rm = TRUE)





dist_df <- crossing(i = 1:nrow(pw_fst@PairwiseFSTmatrix), j = i) %>%
    filter(i > 1, j < i) %>%
    pmap_dfr(function(i, j) {
        si <- which(samp_df$location == colnames(pw_fst@PairwiseFSTmatrix)[i])[1]
        sj <- which(samp_df$location == colnames(pw_fst@PairwiseFSTmatrix)[j])[1]
        .dist <- sqrt((samp_df$lat[si] - samp_df$lat[sj])^2 +
                          (samp_df$lon[si] - samp_df$lon[sj])^2)
        # Convert to km:
        .dist <- .dist / 1000
        .fst <- pw_fst@PairwiseFSTmatrix[i,j]
        .loc1 <- colnames(pw_fst@PairwiseFSTmatrix)[i]
        .loc2 <- colnames(pw_fst@PairwiseFSTmatrix)[j]
        tibble(dist = .dist, fst = .fst, loc1 = .loc1, loc2 = .loc2,
               i1 = i, i2 = j)
    }) %>%
    mutate(n_adults = map2_int(loc1, loc2,
                               function(a, b) {
                                   dd <- filter(samp_df, location %in% c(a,b))
                                   min(dd$n_adults)
                               }))

dist_fst <- pw_fst@PairwiseFSTmatrix
diag(dist_fst) <- 0

dist_km <- matrix(0, nrow(pw_fst@PairwiseFSTmatrix), ncol(pw_fst@PairwiseFSTmatrix))
for (i in 2:nrow(dist_km)) {
    for (j in 1:(i-1)) {
        x <- filter(dist_df, i1 == i, i2 == j)[["dist"]]
        dist_km[i,j] <- dist_km[j,i] <- x
    }
}
# This tests for an effect of distance. Very obviously significant.
km_fst_test <- mantel(dist_km, dist_fst, permutations = 2000)





# int<lower=1> N;                   // # observations
# int<lower=2> P;                   // # populations
# vector[N] fst;                    // pairwise Fst
# vector[N] dist;                   // pairwise distances
# int<lower=1, upper=P> pop1[N];    // population id for first of pair
# int<lower=1, upper=P> pop2[N];    // population id for second of pair


fst_pools <- unique(c(dist_df$loc1, dist_df$loc2))
stan_data <- list(N = nrow(dist_df),
                  P = length(unique(c(dist_df$loc1, dist_df$loc2))),
                  fst = dist_df$fst,
                  dist = dist_df$dist,
                  pop1 = map_int(dist_df$loc1, ~ which(fst_pools == .x)),
                  pop2 = map_int(dist_df$loc2, ~ which(fst_pools == .x)))
# z-score Fst and distances:
stan_data$fst <- (dist_df$fst - mean(dist_df$fst)) / sd(dist_df$fst)
stan_data$dist <- (dist_df$dist - mean(dist_df$dist)) / sd(dist_df$dist)

fst_mod <- stan("fst-dist.stan",
                data = stan_data,
                chains = 4,
                iter = 10e3L)

plot(fst_mod, pars = c("pop_int", "beta"))
plot(fst_mod, plotfun = "hist", pars = c("sig_res", "sig_pop"), ncol = 1)

class(fst_mod@sim[[1]][["alpha[1]"]])


rstan::extract(fst_mod, c("sig_res", "sig_pop")) %>%
    do.call(what = cbind) %>%
    as.data.frame() %>%
    set_names(c("sig_res", "sig_pop")) %>%
    as_tibble() %>%
    pivot_longer(everything()) %>%
    ggplot(aes(value)) +
    geom_freqpoly(aes(color = name), bins = 50)


rstan::extract(fst_mod, c("alpha"))[[1]] %>%
    as.data.frame() %>%
    set_names(c("b0", "b1")) %>%
    as_tibble() %>%
    pivot_longer(everything()) %>%
    ggplot(aes(value)) +
    geom_freqpoly(aes(color = name), bins = 50)



rstan::extract(fst_mod, "pop_int")[[1]] %>%
    as.data.frame() %>%
    set_names(fst_pools) %>%
    as_tibble() %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarize(lo = quantile(value, 0.025),
              med = quantile(value, 0.5),
              hi = quantile(value, 0.975)) %>%
    ggplot(aes(med, name)) +
    geom_vline(xintercept = 0, color = "gray70", linetype = 2) +
    geom_pointrange(aes(xmin = lo, xmax = hi)) +
    theme_classic()

