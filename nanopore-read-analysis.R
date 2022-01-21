
library(tidyverse)
library(progress)
library(Rcpp)
library(RcppProgress)
library(BH)




# ln -s /opt/R/arm64/gfortran/lib/libgfortran.5.dylib /usr/lib/libgfortran.5.dylib

# cp /opt/R/arm64/gfortran/lib/libgfortran.5.dylib /usr/lib/


# Should have boost installed!
Sys.setenv("PKG_LIBS"="-lboost_iostreams")
# Sys.setenv(FLIBS = "-L/usr/local/gfortran/lib/gcc/aarch64-apple-darwin20.2.0/11.0.0 -L/usr/local/gfortran/lib -lgfortran -lm")
# Sys.setenv("CFLAGS" = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
# Sys.setenv("CCFLAGS" = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
# Sys.setenv("CPPFLAGS" = "-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
# Sys.setenv(CPPFLAGS = "-Xclang -fopenmp")
# Sys.setenv(LDFLAGS = "-lomp")
Sys.unsetenv("FLIBS")
Sys.unsetenv("CFLAGS")
Sys.unsetenv("CCFLAGS")
Sys.unsetenv("CPPFLAGS")
Sys.unsetenv("CPPFLAGS")
Sys.unsetenv("LDFLAGS")

Rcpp::sourceCpp("nanopore-read-analysis.cpp", rebuild = TRUE)


# -----------*
# First make sorted version of summary file:
# -----------*

ss <- read_tsv("~/_data/basecalls_guppy-5.0.11--sequencing_summary.txt.gz") %>%
    # Ones that didn't pass filters aren't in FASTQ file anyway:
    filter(passes_filtering)
# Takes ~ 11 min
all_names <- get_read_names(in_fn = "~/_data/basecalls_guppy-5.0.11.fastq.gz",
                            n_reads = nrow(ss))
length(all_names); length(ss$read_id)
head(all_names); head(ss$read_id)
check_names <- all_names == ss$read_id
# Yep, they're different:
mean(!check_names)
which(!check_names)[1]
# They're all in there, though:
find_names <- all_names %in% ss$read_id
mean(find_names)
# Write output:
write_lines(all_names, "~/_data/basecalls_guppy-5.0.11_read_names.txt")
# Make sorted version of summary file:
ss2 <- ss
ss2 <- ss2[match(all_names, ss2$read_id),]
# write_tsv(ss2, "~/_data/basecalls_guppy-5.0.11--sorted_sequencing_summary.txt")


# -----------*
# Now use the sorted summary to filter and make list of ordered line numbers
# This should make extracting reads faster.
# -----------*

ss <- read_tsv("~/_data/basecalls_guppy-5.0.11--sorted_sequencing_summary.txt.gz")

ss2 <- read_tsv("~/_data/sorted_basecalls_guppy-5.0.11--sequencing_summary.txt.gz")

ss2 <- ss2[,-1]
ss
ss2
which(ss$read_id != ss2$read_id)[1]
cbind(head(ss$read_id), head(ss2$read_id))
cbind(ss$read_id[3760], ss2$read_id[3760])




hp <- ss %>%
    ggplot( aes(mean_qscore_template, sequence_length_template)) +
    geom_hex() +
    scale_color_viridis_c()
# hp

# hist(ss$sequence_length_template)
# hist(ss$mean_qscore_template)
# hist(ss$strand_score_template)

lq_filter <- ss$mean_qscore_template > 10 & ss$sequence_length_template > 14.5e3

sum(ss$sequence_length_template[lq_filter]) / 100e6

clean_reads_names <- ss$read_id[lq_filter]
# Read numbers converted to 0-based indices
clean_reads <- which(lq_filter) - 1




# -----------*
# Then do the actual filtering:
# -----------*



# Takes ~17.5 min

filter_fastq(in_fn = "~/_data/basecalls_guppy-5.0.11.fastq.gz",
             out_fn = "~/_data/basecalls_guppy-5.0.11_filtered.fastq.gz",
             read_nums = clean_reads)




all_names_f <- get_read_names(in_fn = "~/_data/basecalls_guppy-5.0.11_filtered.fastq.gz",
                              n_reads = length(clean_reads))

identical(all_names_f, clean_reads_names)



# ## Create connections
# con_in <- gzfile(description="~/_data/basecalls_guppy-5.0.11.fastq.gz", open="r")
# con_out <- gzfile(description="~/_data/basecalls_guppy-5.0.11_filtered.fastq.gz", open="w")
#
# # Total number of reads in FASTQ
# n <- nrow(ss)
#
# pb <- progress_bar$new(
#     format = "  [:bar] :percent eta: :eta",
#     total = n, clear = FALSE, width= 60)
#
# good_read_idx <- 1
# filtered_bp <- 0
# n_filtered_reads <- 0
# for(i in 1:n) {
#     tmp <- scan(file=con_in, nlines=4, quiet=TRUE, what = character(), sep = "\n")
#     if (grepl(paste0("^@", reads[good_read_idx]), tmp[1])) {
#         write(tmp, con_out, append = TRUE, ncolumns = 1)
#         good_read_idx <- good_read_idx + 1
#         filtered_bp <- filtered_bp + nchar(tmp[2])
#         n_filtered_reads <- n_filtered_reads + 1
#     }
#     pb$tick()
# }; print(filtered_bp / 1e9)
#
# close(con_in)
# close(con_out)

