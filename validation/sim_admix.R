#!/usr/bin/env Rscript
# In-silico dilution of a ballad sample with normal cells of the SAME person.
#   alleles: per panel site, germline genotype sampled from its posterior (1kGP prior, het rate h, the fitted tumor
#            state at p = 1 with the theta0 floor); a share w = 2(1-p) / (p c + 2(1-p)) of the site's reads is replaced by
#            normal reads drawn from that genotype (hets at 0.5 +/- theta0, random side); tumor reads are a subsample of
#            the observed reads (no added noise); depth is unchanged.
#   expression: logFC' = log2(p 2^(logFC - base0) + (1 - p)) + base0 (normal cells at the diploid reference level).
# Usage: sim_admix.R <src_dir> <dst_dir> <p> <seed> <h> <theta0> <base0>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
src <- a[1]; dst <- a[2]; p <- as.numeric(a[3]); set.seed(as.integer(a[4]))
h <- as.numeric(a[5]); th0 <- as.numeric(a[6]); base0 <- as.numeric(a[7]); eps <- 0.01
dir.create(file.path(dst, "numbat"), recursive = TRUE, showWarnings = FALSE)

ac <- fread(file.path(src, "allele_counts.tsv.gz"))
sg <- fread(file.path(src, "calls_segments.tsv"))[!is.na(c)]
ac[, chr := as.integer(sub("chr", "", CHROM))]
ac[, `:=`(cc = 2, mm = 1)]
for (i in seq_len(nrow(sg))) ac[chr == sg$CHROM[i] & POS >= sg$start[i] & POS <= sg$end[i], `:=`(cc = sg$c[i], mm = sg$minor[i])]
ac[, th := pmin(sqrt(abs((cc - mm) / cc - 0.5)^2 + th0^2), 0.495)]
f <- ac$AF; d <- ac$DP; x <- ac$ALT_N
lp <- cbind(log((1 - f)^2 + (1 - h) * f * (1 - f)) + dbinom(x, d, eps, log = TRUE),
            log(2 * h * f * (1 - f)) + log(0.5 * (dbinom(x, d, 0.5 + ac$th) + dbinom(x, d, 0.5 - ac$th))),
            log(f^2 + (1 - h) * f * (1 - f)) + dbinom(x, d, 1 - eps, log = TRUE))
pp <- exp(lp - apply(lp, 1, max)); pp <- pp / rowSums(pp)
u <- runif(nrow(pp)); g <- 1L + (u > pp[, 1]) + (u > pp[, 1] + pp[, 2])        # 1 hom-ref, 2 het, 3 hom-alt
w <- 2 * (1 - p) / (p * ac$cc + 2 * (1 - p))
n_norm <- rbinom(nrow(ac), d, w); n_tum <- d - n_norm
alt_t <- ifelse(n_tum > 0, rhyper(nrow(ac), x, d - x, n_tum), 0L)
af_n <- c(eps, NA, 1 - eps)[g]
af_n[g == 2] <- 0.5 + th0 * sample(c(-1, 1), sum(g == 2), replace = TRUE)
alt_n <- rbinom(nrow(ac), n_norm, af_n)
ac[, ALT_N := alt_t + alt_n][, REF_N := DP - ALT_N]
fwrite(ac[, .(CHROM, POS, REF, ALT, AF, REF_N, ALT_N, DP, OTH)], file.path(dst, "allele_counts.tsv.gz"), sep = "\t")

bk <- fread(file.path(src, "numbat", "bulk_joint.tsv.gz"))
bk[!is.na(logFC), logFC := log2(p * 2^(logFC - base0) + (1 - p)) + base0]
fwrite(bk, file.path(dst, "numbat", "bulk_joint.tsv.gz"), sep = "\t")
file.copy(file.path(src, c("karyotype.tsv", "summary.tsv")), dst, overwrite = TRUE)
cat(sprintf("p = %.2f: %d sites, genotypes hom-ref/het/hom-alt %s, mean normal-read share %.3f\n",
            p, nrow(ac), paste(tabulate(g, 3), collapse = "/"), mean(w)))
