# Genotype-free allelic-imbalance model shared by plot_genome.R and call_arms.R.
#
# Every panel SNP with reads (any depth) is a mixture of hom-ref / het / hom-alt with 1kGP prior
# (1-f)^2 + (1-h)f(1-f) / 2hf(1-f) / f^2 + (1-h)f(1-f): h < 1 lets RNA show fewer hets than 1kGP predicts
# (mono-allelic expression, ancestry, homozygous runs) and is fitted per sample. Het reads are
# 0.5 Bin(d, 0.5 + theta) + 0.5 Bin(d, 0.5 - theta): no phasing (RNA phase is unreliable beyond ~50 kb) and no genotype
# calling (calling first drops low-depth hets that show one allele and hides LOH). Homs: Bin(d, eps) / Bin(d, 1 - eps).
# Bins close when the expected number of informative hets, sum 2f(1-f)(1 - 0.5^(d-1)), reaches bin_snps. Each bin keeps
# its log-likelihood curve over theta_grid, so segments are fitted by adding curves rather than averaging estimates.
# Tested alternatives that did not help on RNA (see experimental/): GLIMPSE2 genotype prior, phase within 100 kb blocks,
# beta-binomial overdispersion (not identifiable without genotypes, and still degenerate with GLIMPSE genotypes).

suppressPackageStartupMessages(library(data.table))

theta_grid <- seq(0, 0.49, by = 0.005)
baf_eps <- 0.01
hg38_len <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
              133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616,
              64444167, 46709983, 50818468)
# centromere midpoints (UCSC cytoband acen, hg38); 13, 14, 15, 21, 22 are acrocentric (no p-arm data)
hg38_cen <- c(123400000, 93900000, 90900000, 50000000, 48750000, 60550000, 60100000, 45200000, 43850000, 39800000,
              53400000, 35500000, 17700000, 17150000, 19000000, 36850000, 25050000, 18450000, 26150000, 28050000,
              11950000, 15550000)
acrocentric <- c(13, 14, 15, 21, 22)

lse <- function(...) { m <- pmax(...); m + log(Reduce(`+`, lapply(list(...), function(v) exp(v - m)))) }

read_sites <- function(sd) {
  # panel SNPs with reads, chr1-22, MHC masked as numbat does
  fread(file.path(sd, "allele_counts.tsv.gz"))[, CHROM := as.integer(sub("chr", "", CHROM))][
    CHROM %in% 1:22 & !(CHROM == 6 & POS > 28.5e6 & POS < 33.5e6) & DP >= 1][order(CHROM, POS)]
}

make_bins <- function(s, bin_snps) {
  # adds s$bin (per chromosome) and returns bins; drops a chromosome-end remainder with < half the target
  s[, w := 2 * AF * (1 - AF) * (1 - 0.5^(DP - 1))][, bin := floor(cumsum(w) / bin_snps), by = CHROM]
  b <- s[, .(start = min(POS), end = max(POS), POS = as.numeric(median(POS)), n_sites = .N, exp_hets = sum(w)),
         by = .(CHROM, bin)][exp_hets >= bin_snps / 2]
  s[b, on = .(CHROM, bin), keep := TRUE]
  s[is.na(keep), bin := NA][, keep := NULL]
  b[order(CHROM, bin)]
}

site_ll <- function(a, d, f, h) {
  # site x theta log-likelihood matrix
  hom <- lse(log((1 - f)^2 + (1 - h) * f * (1 - f)) + dbinom(a, d, baf_eps, log = TRUE),
             log(f^2 + (1 - h) * f * (1 - f)) + dbinom(a, d, 1 - baf_eps, log = TRUE))
  lh <- log(2 * h * f * (1 - f)) + log(0.5)
  vapply(theta_grid, function(t)
    lse(hom, lh + lse(dbinom(a, d, 0.5 + t, log = TRUE), dbinom(a, d, 0.5 - t, log = TRUE))), numeric(length(a)))
}

bin_curves <- function(s, b, h) {
  # bins x theta log-likelihood curves, row order = b
  x <- s[!is.na(bin)]
  m <- do.call(rbind, lapply(split(x, x$CHROM), function(y) {
    r <- rowsum(site_ll(y$ALT_N, y$DP, y$AF, h), y$bin)
    rownames(r) <- paste(y$CHROM[1], rownames(r)); r
  }))
  m[paste(b$CHROM, b$bin), , drop = FALSE]
}

fit_h <- function(s, b, n_fit = 300, h_grid = seq(0.2, 1, by = 0.05)) {
  # sample-wide het rate: maximise the sum over a sample of bins of each bin's best theta
  set.seed(1)
  fb <- b[sample(.N, min(n_fit, .N))]
  x <- s[fb, on = .(CHROM, bin), nomatch = NULL]
  ll <- vapply(h_grid, function(h) {
    r <- rowsum(site_ll(x$ALT_N, x$DP, x$AF, h), paste(x$CHROM, x$bin))
    sum(apply(r, 1, max))
  }, 0)
  h_grid[which.max(ll)]
}

theta_fit <- function(curve) {
  # MLE on the grid and a 1-SE support interval (log-likelihood within 0.5 of the maximum)
  i <- which.max(curve); ok <- which(curve >= curve[i] - 0.5)
  list(theta = theta_grid[i], lo = theta_grid[min(ok)], hi = theta_grid[max(ok)],
       se = max((theta_grid[max(ok)] - theta_grid[min(ok)]) / 2, 0.005))
}

allele_bins <- function(sd, bin_snps = 20) {
  # full model for one sample: sites, bins (with theta), curves, h
  s <- read_sites(sd)
  b <- make_bins(s, bin_snps)
  h <- fit_h(s, b)
  cv <- bin_curves(s, b, h)
  f <- lapply(seq_len(nrow(cv)), function(i) theta_fit(cv[i, ]))
  b[, `:=`(theta = vapply(f, `[[`, 0, "theta"), theta_se = vapply(f, `[[`, 0, "se"))]
  list(sites = s, bins = b, curves = cv, h = h)
}
