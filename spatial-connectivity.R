
library(tidyverse)
library(lubridate)
library(readxl)
library(rgdal)
library(broom)
library(patchwork)
library(viridis)
library(poolfstat)

source("helpers.R")



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
    filter(year(date) == 2019) %>%
    select(biotech_id, contam_pilot, n_adults, date, lake, site, location, lat, lon) %>%
    to_utm()


pools <- read_lines("~/_data/snape/space_sync_masked.names.gz") %>%
    map_chr(~ samp_df$location[samp_df$biotech_id == .x])
n_adults <- map_int(pools, ~ samp_df$n_adults[samp_df$location == .x][1])

pool_dat <- popsync2pooldata(sync.file = "~/_data/snape/space_sync_masked_noblanks.sync.gz",
                             poolsizes = n_adults * 2,
                             poolnames = pools,
                             min.maf = 0,
                             nthreads = 4)

# Identity method seems to work better for dealing with low number of
# individuals in KS pool
pw_fst <- compute.pairwiseFST(pool_dat, method = "Identity",
                               nsnp.per.bjack.block = 0)
                               # nsnp.per.bjack.block = pool_dat@nsnp %/% 500)
plot(pw_fst, cex = 0.5)

fst_mat <- pw_fst@PairwiseFSTmatrix
fst_mat


one_loc_fct <- function(lake, .loc1, .loc2) {
    factor(as.integer(.loc1 == lake | .loc2 == lake),
           levels = 0:1, labels = c("others", lake))
}

dist_df <- crossing(i = 1:nrow(fst_mat), j = i) %>%
    filter(i > 1, j < i) %>%
    pmap_dfr(function(i, j) {
        si <- which(samp_df$location == colnames(fst_mat)[i])[1]
        sj <- which(samp_df$location == colnames(fst_mat)[j])[1]
        .dist <- sqrt((samp_df$lat[si] - samp_df$lat[sj])^2 +
                                    (samp_df$lon[si] - samp_df$lon[sj])^2)
        # Convert to km:
        .dist <- .dist / 1000
        .fst <- fst_mat[i,j]
        .loc1 <- colnames(fst_mat)[i]
        .loc2 <- colnames(fst_mat)[j]
        tibble(dist = .dist, fst = .fst, loc1 = .loc1, loc2 = .loc2)
    }) %>%
    mutate(tjor = one_loc_fct("Tjörnin", loc1, loc2),
           n_adults = map2_int(loc1, loc2,
                               function(a, b) {
                                   dd <- filter(samp_df, location %in% c(a,b))
                                   min(dd$n_adults)
                               }))

other_dist_df <- dist_df[dist_df$tjor == "others",] %>% select(-tjor)

dist_p_list <- c("all", "Elliðavatn", "Ólafsfjarðarvatn", "Blikalón",
                 "Miklavatn-W", "Skjálftavatn", "Þernuvatn", "Tjörnin") %>%
    set_names() %>%
    map(function(pp) {
        if (pp == "all" || pp == "Tjörnin") {
            dd <- other_dist_df %>%
                mutate(site = factor(rep(1, n())))
            if (pp == "all") pp <- ""
        } else {
            dd <- other_dist_df %>%
                mutate(site = one_loc_fct(pp, loc1, loc2)) %>%
                arrange(site)
        }
        p <- dd %>%
            ggplot(aes(dist, fst, color = site)) +
            ggtitle(pp) +
            geom_point(size = 2) +
            stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
            xlab("Distance (km)") +
            ylab(expression(F[ST])) +
            scale_color_manual(NULL, values = c("gold2", "dodgerblue3"),
                               guide = "none") +
            theme_minimal() +
            NULL
        if (pp == "Tjörnin") {
            tdd <- dist_df[dist_df$tjor == "Tjörnin",]
            p <- p +
                geom_point(data = tdd, size = 2, color = "dodgerblue3") +
                stat_smooth(data = tdd, method = "lm", formula = y ~ x,
                            se = FALSE, color = "dodgerblue3") +
                ggtitle("Tjörnin")
        }
        return(p)
    })
wrap_plots(dist_p_list, nrow = 2)

for (x in c("all", "Elliðavatn", "Ólafsfjarðarvatn", "Blikalón", "Tjörnin")) {
    ggsave(sprintf("_plots/dist_fst_%s.pdf", x), dist_p_list[[x]],
           width = 3, height = 3)
}







dist_fst_p

dist_fst_p +
    geom_point(data = dist_df[dist_df$tjor == "Tjörnin",],
               color = "black", fill = "gray90", shape = 21, size = 4)



lm(fst ~ dist, dist_df[dist_df$tjor == "Tjörnin",]) %>% summary()
lm(fst ~ dist + elli, other_dist_df) %>% summary()

library(lme4)
z <- lmer(fst ~ dist + (1 | loc1) + (1 | loc2),
          other_dist_df %>%
              mutate(dist = (dist - mean(dist)) / sd(dist)))

z %>% summary()

get_ranef <- function(x) {
    rn <- rownames(ranef(x)[["loc1"]])[rownames(ranef(x)[["loc1"]]) %in%
                                           rownames(ranef(x)[["loc2"]])]
    o <- ranef(x)[["loc1"]][rn,] + ranef(x)[["loc2"]][rn,]
    names(o) <- rn
    return(o)
}

bz1 <- bootMer(z, get_ranef, nsim = 2000, use.u = TRUE)
bz2 <- bootMer(z, function(x) fixef(x)[["dist"]], nsim = 2000, use.u = FALSE)
bzm1 <- apply(bz1$t, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% t()
quantile(bz2$t, probs = c(0.025, 0.5, 0.975))

# bzm1[(bzm1[,1] < 0 & bzm1[,3] < 0) | (bzm1[,1] > 0 & bzm1[,3] > 0),] %>%
bzm1 %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c("loc", "lo", "med", "hi")) %>%
    ggplot(aes(med, loc)) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray70", size = 0.75) +
    geom_pointrange(aes(xmin = lo, xmax = hi), size = 1, fatten = 2) +
    theme_minimal()







# sim.graph <- map(pools, ~ c(.x, "A", "")) %>%
sim.graph <- map(LETTERS[1:7], ~ c(.x, "R", "")) %>%
    do.call(what = rbind)
sim.graph.params <- generate.graph.params(sim.graph)
sim.graph.params %>% plot()








spat_af <- read_tsv("~/_data/snape/space_snape_masked_noblanks.sync.gz",
                    col_names = c("scaff", "pos", "ref", pools),
                    col_types = paste(c("cic", rep("c", length(pools))),
                                      collapse = ""))

# Kapun et al. (2021) only used SNPs > 500 bp away from each other.
spat_af <- spat_af %>%
    filter(min_dist_filter(scaff, pos, 500))

for (p in pools) {
    spat_af[[p]] <- spat_af[[p]] %>%
        str_split(":") %>%
        map_chr(~ .x[[7]]) %>%
        as.numeric()
}

af_mat <- spat_af %>%
    # .[,-which(colnames(spat_af) %in% samp_df$biotech_id[samp_df$contam_pilot])] %>%
    .[,-1:-3] %>%
    as.matrix() %>%
    t()


af_pc <- prcomp(af_mat)

af_pc %>% summary()
af_pc %>% str()

# Cluster analysis:
wss <- (nrow(af_pc$x)-1)*sum(apply(af_pc$x,2,var))
for (i in 2:(nrow(af_pc$x)-1)) {
    wss[i] <- sum(kmeans(af_pc$x, centers=i, nstart=25, iter.max=1000)$withinss)
}
plot(1:length(wss), wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
k <- kmeans(af_pc$x, 4, nstart=25, iter.max=1000)


af_pc_df <- af_pc$x %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    as_tibble() %>%
    mutate(location = map_chr(sample, ~ samp_df$location[samp_df$biotech_id == .x][1]),
           lat = map_dbl(sample, ~ samp_df$lat[samp_df$biotech_id == .x][1]),
           lon = map_dbl(sample, ~ samp_df$lon[samp_df$biotech_id == .x][1]),
           clust = map_int(sample, ~ k$clust[[.x]]) %>%
               factor())

library(ggrepel)

af_pc_df %>%
    # mutate(location = factor(location, levels = c(
    #     "Áshildarholtsvatn", "Blikalón", "Grafarvatn", "Elliðavatn",
    #     "Hrisatjorn", "Leirhafnarvatn", "Lýsuvatn", "Miklavatn",
    #     "Miklavatn-W", "Fellshóll", "Kálfaströnd", "Syðri Neslönd",
    #     "Ólafsfjarðarvatn", "Þernuvatn", "Skjálftavatn", "Svartárvatn",
    #     "Tjörnin", "Víkingavatn"))) %>%
    ggplot(aes(PC1, PC2, color = clust)) +
    geom_point() +
    # geom_text(aes(label = location, nudge_x = nx, hjust = hj),
    #           size = 6 / 2.83465) +
    geom_text_repel(aes(label = location)) +
    # scale_color_manual(values = rep(viridis(4, end = 0.8), 5)) +
    scale_color_viridis_d(begin = 0.2) +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(NA, max(af_pc_df$PC1) * 1.1)) +
    NULL



# ============================================================================*
# ============================================================================*
# ============================================================================*
# ============================================================================*





myvatn_df <- readOGR(dsn = paste0("~/Box Sync/midgenomics/location_data/",
                                  "shapefiles/myvatn"),
                     layer = "Myvatn_WSGUTM28") %>%
    tidy() %>%
    rename(lon = long)
# Iceland outline is from GADM data (version 3.6; https://gadm.org/)
iceland_df <- readOGR(dsn = paste0("~/Box Sync/midgenomics/location_data/",
                                   "shapefiles/iceland"),
                      layer = "gadm36_ISL_0") %>%
    tidy() %>%
    rename(lon = long) %>%
    to_utm() %>%
    # -----------`
    # Shown below are two ways to filter this dataset, to avoid
    # plotting islands far from shore:
    # -----------`
    # 1. You can filter out islands that are very far from shore:
    # filter(!piece %in% c(90, 133, 143, 157, 215, 244, 257, 258, 260, 262))
    # 2. Filter for just the mainland:
    filter(piece == 1)

# lakes_iceland_p <-
iceland_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group), color = "black", size = 0.1,
                 fill = "gray80") +
    geom_point(data = af_pc_df, aes(color = clust), size = 6) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 0.75, fill = NA, color = "black") +
    coord_equal() +
    # scale_color_viridis_c(end = 0.9) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL


lakes_myvatn_p <- myvatn_df %>%
    ggplot(aes(x = lon, y = lat)) +
    geom_rect(data = tibble(x = 403118.1 - 100, xe = 411932.0 + 100,
                            y = 7271491 - 100, ye = 7282704 + 100),
              aes(xmin = x, xmax = xe, ymin = y, ymax = ye),
              inherit.aes = FALSE, size = 0.75, fill = "white", color = "black") +
    geom_polygon(aes(group = group, fill = hole), color = "black", size = 0.1) +
    geom_point(data = samp_df %>%
                   filter(biotech_id %in% pools, lake == "Mývatn"),
               size = 6, color = "dodgerblue3") +
    geom_segment(data = tibble(lone = 405e3, lon = lone - 1e3, lat = 7280e3),
                 aes(xend = lone, yend = lat), size = 1) +
    geom_text(data = tibble(lat = 7280e3 + 200,  lon = 405e3-500),
              label = "1 km", vjust = 0, hjust = 0.5, size = 10 / 2.83465) +
    scale_fill_manual(values = c("lightblue1", "white"), guide = "none") +
    coord_equal(xlim = c(403118.1 - 100, 411932.0 + 100),
                ylim = c(7271491 - 100, 7282704 + 100)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL

ggsave("~/Desktop/pilots-iceland.pdf", lakes_iceland_p, width = 6, height = 4)


ggsave("~/Desktop/pilots-myvatn.pdf", lakes_myvatn_p, width = 2.4, height = 3)








# ============================================================================*
# ============================================================================*

# Convert to MakeTree format ----

# ============================================================================*
# ============================================================================*


library(tidyverse)
library(ape)

source("helpers.R")

# This sets plotting device on my computer:
if (file.exists(".Rprofile")) source(".Rprofile")

pools <- read_lines("~/_data/snape/space_sync_masked.names.gz")

# Allele frequencies
# Takes ~ 9 sec
freq_df <- read_tsv( "~/_data/snape/space_snape_masked_noblanks.sync.gz",
                     col_names = c("scaff", "pos", "ref", pools),
                     col_types = paste(c("cic", rep("c", length(pools))),
                                       collapse = ""),
                     progress = FALSE) %>%
    mutate(across(all_of(pools), get_snape_af)) %>%
    correct_biallelic_freqs("~/_data/snape/space_sync_masked_noblanks.sync.gz")

# # Allele counts
# # Takes ~9 sec
# count_df <- read_tsv( "~/_data/snape/space_sync_masked_noblanks.sync.gz",
#                       col_names = c("scaff", "pos", "ref", pools),
#                       col_types = paste(c("cic", rep("c", length(pools))),
#                                         collapse = ""),
#                       progress = FALSE) %>%
#     mutate(across(all_of(pools), split_sync_strings)) %>%
#     identity()


#' The allele frequencies are in the format A:T:C:G:N:del,
#' i.e: count of bases 'A', count of bases 'T',... and deletion count
#' in the end













allele_strings <- maketree_format(allele_counts) %>%
    set_names(pools) %>%
    as_tibble() %>%
    filter(if_all(.fns = ~ .x != ""))


write_delim(allele_strings, "~/_data/snape/space_sync_masked_biallelic.snp.gz", delim = " ")


source("plotting_funcs.R")


plot_tree("~/_data/snape/admixture/m0/m0_s1193447534")
plot_resid("~/_data/snape/admixture/m0/space_m0", gzfile("~/_data/snape/space_sync_masked.names.gz"))





"(((Vik-19_S41:0.00199578,(Mik-19_S15:0.00191953,(((MyBR-19_S17:0.00249207,MyKS-19-B_S18:0.00233249):4.67111e-05,MySN-19_S19:0.00235689):0.000133806,Sva-19_S39:0.00240741):0.000977135):0.000457346):0.000283039,Skj-19_S22:0.0018737):0.00164882,(((Lei-19_S13:0.00255718,Blik-19_S6:0.00914425):0.00126967,Pernuvatn_S21:0.0101247):0.00182329,((MikW-19_S16:0.00438326,(Ola-19_S20:0.00754027,Tjornin_S40:0.0164645):0.000722681):0.00312537,(Hrisatjorn_S11:0.00388747,((Blo-19_S7:0.00321111,Ash-19_S5:0.00371718):0.000164562,(Lys-19_S14:0.00738896,Ellid-19_S8:0.00732375):0.00296746):0.00265435):0.00049805):0.00160464):0.00164882);" %>%
    {ape::read.tree(text = .)} %>%
    plot()


trees <- c("1193447534", "326161346", "619752663", "740815883", "800249271") %>%
    sprintf(fmt = "~/_data/snape/admixture/m0/m0_s%s.treeout.gz") %>%
    map(read.tree) %>%
    do.call(what = c)


par(mfrow = c(2, 3))
for (i in 1:length(trees)) {
    plot(trees[[i]])
}
plot.new()  # just a blank to avoid overflow

trees2 <- "1853570151" %>%
    sprintf(fmt = "~/_data/snape/admixture/m0/m0_s%s.treeout.gz") %>%
    read.tree()

par(mfrow = c(1, 1))
trees2 %>% plot()


# LogLik = -66.374812236
best <- c("1367371382", "592787828", "210721078", "1967952086", "60956383", "1248603299", "1294051240", "838607120", "2040667896")

# LogLik = -106.460367458
worst <- "808069704"


trees3 <- best %>%
    sprintf(fmt = "~/_data/snape/admixture/m0/m0_s%s.treeout.gz") %>%
    map(read.tree) %>%
    do.call(what = c)

par(mfrow = c(3, 3))
for (i in 1:9) plot(trees3[[i]])
# for (i in 1:9) heatmap(trees3[[i]])


par(mfrow = c(1, 1))
"(((Vik-19_S41:0.00199578,(Mik-19_S15:0.00191953,((MySN-19_S19:0.00235689,(MyBR-19_S17:0.00249207,MyKS-19-B_S18:0.00233249):4.67111e-05):0.000133806,Sva-19_S39:0.00240741):0.000977135):0.000457346):0.000283039,Skj-19_S22:0.0018737):0.00176238,(((Lei-19_S13:0.00276826,Pernuvatn_S21:0.00922479):0.00182379,((MikW-19_S16:0.00438326,(Ola-19_S20:0.00754027,Tjornin_S40:0.0164645):0.000722681):0.00131691,Blik-19_S6:0.0100658):0.00141663):0.00194216,(Hrisatjorn_S11:0.00388747,(Blo-19_S7:0.00341019,(Ash-19_S5:0.00369561,(Ellid-19_S8:0.00732375,Lys-19_S14:0.00738896):0.00277825):0.000294186):0.00256622):0.00101277):0.00176238);" %>%
    {read.tree(text = .)} %>%
    plot()
