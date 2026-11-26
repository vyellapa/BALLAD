#!/usr/bin/env Rscript
# Experimental: does a GLIMPSE2 genotype prior make RNA BAF bins better than the 1kGP (HWE) prior?
# Needs a finished sample directory plus <sample_dir>/glimpse/<sample>.chr<N>.bcf from experimental/glimpse_rna.sh.
#
# Unphased BAF bins (theta, drawn at 0.5 +/- theta), same model as R/plot_genome.R except the genotype prior:
#   HWE      P(G) = (1-f)^2, 2f(1-f), f^2 from 1kGP AF, het prior scaled by a fitted sample-wide h
#   GLIMPSE  P(G) from GLIMPSE2's genotype posterior GP with this site's own read likelihood divided back out
#            (GP_k / L_k, L under a balanced het), so the reads are not counted twice and the prior is LD only
#   GLIMPSE + beta-binomial  same prior, het reads beta-binomial (= each allele's count negative binomial with shared
#            dispersion) with concentration g fitted on a sample of bins: allele-specific expression / mapping bias
#            make RNA allele ratios overdispersed, so deep sites stop dominating. Without a genotype prior this was not
#            identifiable (R/plot_genome.R --overdispersion); with GLIMPSE it should be.
# Locally phased + genotype-weighted: same sites and genotype prior. RNA hets cluster in expressed genes with large
#   gaps, so phase is trusted only within blocks (GLIMPSE-phased hets with no gap > --block_kb): inside a block the hap1
#   read count ~ Bin(d, 0.5 + s theta) with one orientation s = +/-1 per block, mixed 50:50 per block, so 1-2 read hets
#   in the same gene add information while switches between genes cost nothing. Unphased sites are their own block.
#   Drawn at 0.5 +/- theta.
# Phase reliability: on chromosomes numbat calls gained, neighbouring hets with depth >= 20 should carry the excess on
#   the same haplotype; concordance by distance (0.5 = random) for Eagle and GLIMPSE -> <sample>_phase_concordance.tsv
# Phased BAF bins (beta-binomial-free ratio sum(hap1 reads) / sum(depth)):
#   Eagle    numbat's pAD (Eagle2 phase, hets called by numbat's rule)
#   GLIMPSE  GLIMPSE2 phased GT at sites it calls het with GP_het >= --min_gp
# Bins: same expected-het bins for both unphased tracks; --bin_snps hets for the phased tracks.
# Output: <sample_dir>/<sample>_baf_glimpse.png and <sample>_baf_glimpse_bins.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(glue)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character"),
  make_option("--sample_dir", type = "character"),
  make_option("--bin_snps", type = "integer", default = 20),
  make_option("--min_gp", type = "double", default = 0.9),
  make_option("--prior_floor", type = "double", default = 1e-3),
  make_option("--block_kb", type = "double", default = 100)
)))
sd <- opt$sample_dir
eps <- 0.01

hg38 <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
          135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
          46709983, 50818468)
offset <- c(0, cumsum(as.numeric(hg38)))[1:22]

ac <- fread(file.path(sd, "allele_counts.tsv.gz"))[, CHROM := as.integer(sub("chr", "", CHROM))][
  CHROM %in% 1:22 & !(CHROM == 6 & POS > 28.5e6 & POS < 33.5e6) & DP >= 1]
gl <- rbindlist(lapply(1:22, function(c) {
  f <- file.path(sd, "glimpse", glue("{opt$sample}.chr{c}.bcf"))
  fread(cmd = glue("bcftools view -i 'FMT/DS>=0' {f} -Ou | bcftools query -f '%POS\\t%REF\\t%ALT\\t[%GT\\t%GP]\\n'"),
        header = FALSE, col.names = c("POS", "REF", "ALT", "GT", "GP"))[, CHROM := c]
}))
gl[, c("gp0", "gp1", "gp2") := tstrsplit(GP, ",", type.convert = TRUE)][, GP := NULL]
d <- merge(ac, gl, by = c("CHROM", "POS", "REF", "ALT"))[order(CHROM, POS)]
message(sprintf("sites with reads: %d, matched to GLIMPSE output: %d", nrow(ac), nrow(d)))

lse <- function(...) { m <- pmax(...); m + log(Reduce(`+`, lapply(list(...), function(v) exp(v - m)))) }
# LD-only prior: divide the read likelihood (balanced het) out of the posterior, floor, renormalise
l0 <- dbinom(d$ALT_N, d$DP, eps, log = TRUE); l1 <- dbinom(d$ALT_N, d$DP, 0.5, log = TRUE)
l2 <- dbinom(d$ALT_N, d$DP, 1 - eps, log = TRUE)
lp <- cbind(log(pmax(d$gp0, 1e-7)) - l0, log(pmax(d$gp1, 1e-7)) - l1, log(pmax(d$gp2, 1e-7)) - l2)
lp <- lp - apply(lp, 1, max)
pr <- exp(lp); pr <- pr / rowSums(pr); pr <- pmax(pr, opt$prior_floor); pr <- pr / rowSums(pr)
d[, `:=`(p0 = pr[, 1], p1 = pr[, 2], p2 = pr[, 3])]

theta_grid <- seq(0, 0.49, by = 0.005)
h_grid <- seq(0.2, 1, by = 0.05)
dbb <- function(x, n, p, g) lchoose(n, x) + lbeta(x + p * g, n - x + (1 - p) * g) - lbeta(p * g, (1 - p) * g)
ll_theta <- function(a, n, P0, P1, P2, g = Inf) {
  # bin log-likelihood over theta_grid for given genotype priors (vectors over sites); g = Inf: binomial hets
  r <- dbinom(a, n, eps, log = TRUE); s <- dbinom(a, n, 1 - eps, log = TRUE)
  dh <- if (is.infinite(g)) function(p) dbinom(a, n, p, log = TRUE) else function(p) dbb(a, n, p, g)
  vapply(theta_grid, function(t)
    sum(lse(log(P0) + r, log(P2) + s, log(P1) + log(0.5) + lse(dh(0.5 + t), dh(0.5 - t)))), 0)
}
d[, w := 2 * AF * (1 - AF) * (1 - 0.5^(DP - 1))][, bin := floor(cumsum(w) / opt$bin_snps), by = CHROM]
keep <- d[, .(ok = sum(w) >= opt$bin_snps / 2), by = .(CHROM, bin)][ok == TRUE]
d <- d[keep, on = .(CHROM, bin)]

# HWE prior with fitted h (as plot_genome.R)
hwe <- function(h) d[, .(P0 = (1 - AF)^2 + (1 - h) * AF * (1 - AF), P1 = 2 * h * AF * (1 - AF),
                         P2 = AF^2 + (1 - h) * AF * (1 - AF))]
set.seed(1); fit <- unique(d[, .(CHROM, bin)])[sample(.N, min(300, .N))]
h_ll <- vapply(h_grid, function(h) {
  x <- cbind(d, hwe(h))[fit, on = .(CHROM, bin)]
  x[, .(m = max(ll_theta(ALT_N, DP, P0, P1, P2))), by = .(CHROM, bin)][, sum(m)]
}, 0)
h_hat <- h_grid[which.max(h_ll)]
d <- cbind(d, hwe(h_hat))
g_grid <- c(5, 10, 20, 40, 80, 160, 320, 640, Inf)
g_ll <- vapply(g_grid, function(g) d[fit, on = .(CHROM, bin)][
  , .(m = max(ll_theta(ALT_N, DP, p0, p1, p2, g))), by = .(CHROM, bin)][, sum(m)], 0)
g_hat <- g_grid[which.max(g_ll)]
message("beta-binomial fit (GLIMPSE prior), log-lik by g: ", paste(sprintf("%g:%.0f", g_grid, g_ll), collapse = " "))
# locally phased, genotype-weighted: hap1 = the haplotype carrying ALT when GT is 0|1 (as in the pBAF GLIMPSE track)
d[, `:=`(phased = GT %in% c("0|1", "1|0"), x1 = fifelse(GT == "1|0", REF_N, ALT_N))]
d[, block := {
  gap <- c(TRUE, diff(POS) > opt$block_kb * 1e3)
  ph_run <- cumsum(gap | !phased | !shift(phased, fill = FALSE))  # a new block at each gap or unphased site
  ph_run
}, by = CHROM]
ll_blocks <- function(a, x1, n, ph, blk, P0, P1, P2) {
  r <- dbinom(a, n, eps, log = TRUE); s <- dbinom(a, n, 1 - eps, log = TRUE); hom <- lse(log(P0) + r, log(P2) + s)
  vapply(theta_grid, function(t) {
    up <- ifelse(ph, dbinom(x1, n, 0.5 + t, log = TRUE), dbinom(a, n, 0.5 + t, log = TRUE))
    dn <- ifelse(ph, dbinom(x1, n, 0.5 - t, log = TRUE), dbinom(a, n, 0.5 - t, log = TRUE))
    bu <- rowsum(lse(hom, log(P1) + up), blk)[, 1]; bd <- rowsum(lse(hom, log(P1) + dn), blk)[, 1]
    sum(log(0.5) + lse(bu, bd))
  }, 0)
}
bins <- d[, .(POS = as.numeric(median(POS)), n_sites = .N, n_blocks = uniqueN(block),
              theta_blocks = theta_grid[which.max(ll_blocks(ALT_N, x1, DP, phased, block, p0, p1, p2))],
              theta_hwe = theta_grid[which.max(ll_theta(ALT_N, DP, P0, P1, P2))],
              theta_glimpse = theta_grid[which.max(ll_theta(ALT_N, DP, p0, p1, p2))],
              theta_glimpse_bb = theta_grid[which.max(ll_theta(ALT_N, DP, p0, p1, p2, g_hat))]), by = .(CHROM, bin)]

# phased bins
bulk <- fread(file.path(sd, "numbat", "bulk_joint.tsv.gz"))
eagle <- bulk[!is.na(pAD) & DP >= 1, .(CHROM, POS, h1 = pAD, DP)][order(CHROM, POS)]
gh <- d[gp1 >= opt$min_gp & GT %in% c("0|1", "1|0"), .(CHROM, POS, h1 = fifelse(GT == "0|1", ALT_N, REF_N), DP)]
pbin <- function(x, src) x[, bin := (seq_len(.N) - 1) %/% opt$bin_snps, by = CHROM][
  , if (.N >= opt$bin_snps / 2) .(POS = as.numeric(median(POS)), y = sum(h1) / sum(DP)), by = .(CHROM, bin)][, src := src]
pb <- rbind(pbin(eagle, "pBAF Eagle"), pbin(gh, "pBAF GLIMPSE"))
message(sprintf("h = %.2f; phased hets: Eagle %d, GLIMPSE (GP_het >= %.2f) %d", h_hat, nrow(eagle), opt$min_gp, nrow(gh)))

# numbat joint state per bin for colouring
segs <- fread(file.path(sd, "numbat", "segments_joint.tsv"))
segs[, state := fifelse(is.na(cnv_state_post) | cnv_state_post == "", "neu", cnv_state_post)]
state_at <- function(x) {
  s <- segs[, .(CHROM, start = seg_start, end = seg_end, state)]; setkey(s, CHROM, start, end)
  x <- copy(x)[, `:=`(start = POS, end = POS)]
  foverlaps(x, s, by.x = c("CHROM", "start", "end"), mult = "first")[, state := fifelse(is.na(state), "neu", state)]$state
}
bins[, state := state_at(bins)]; pb[, state := state_at(pb)]
fwrite(bins, file.path(sd, glue("{opt$sample}_baf_glimpse_bins.tsv")), sep = "\t")

# spread of theta in neutral bins (noise) vs gained bins (signal: 1/6 expected for a clonal single-copy gain)
summ <- melt(bins, id.vars = c("CHROM", "bin", "state"),
             measure.vars = c("theta_hwe", "theta_glimpse", "theta_glimpse_bb", "theta_blocks"))[
  , .(bins = .N, median_theta = round(median(value), 3), iqr = round(IQR(value), 3),
      at_edge = sum(value >= 0.45)), by = .(prior = variable, state)][order(state, prior)]
print(summ)

# phase concordance on gained chromosomes: does the excess sit on the same haplotype at neighbouring deep hets?
gained <- segs[state %in% c("amp"), unique(CHROM)]
conc <- function(x, src) {
  x <- x[CHROM %in% gained & DP >= 20][order(CHROM, POS)][, major := sign(h1 / DP - 0.5)][major != 0]
  x[, `:=`(dist = shift(POS, -1) - POS, same = major == shift(major, -1)), by = CHROM][!is.na(dist)][
    , .(pairs = .N, concordance = round(mean(same), 3)),
    by = .(distance = cut(dist, c(0, 1e4, 5e4, 2e5, 1e6, Inf), labels = c("<10kb", "10-50kb", "50-200kb", "200kb-1Mb", ">1Mb")))][
    order(distance)][, phasing := src]
}
pc <- rbind(conc(eagle, "Eagle"), conc(gh, "GLIMPSE"))
print(dcast(pc, distance ~ phasing, value.var = c("concordance", "pairs")))
fwrite(pc, file.path(sd, glue("{opt$sample}_phase_concordance.tsv")), sep = "\t")

panels <- c("BAF, HWE prior", "BAF, GLIMPSE prior", glue("BAF, GLIMPSE + beta-binomial (g = {g_hat})"),
            "pBAF Eagle", "pBAF GLIMPSE", glue("BAF, GLIMPSE prior, phased within {opt$block_kb} kb blocks"))
pts <- rbind(bins[, .(CHROM, POS, y = 0.5 + theta_hwe, state, panel = panels[1])],
             bins[, .(CHROM, POS, y = 0.5 - theta_hwe, state, panel = panels[1])],
             bins[, .(CHROM, POS, y = 0.5 + theta_glimpse, state, panel = panels[2])],
             bins[, .(CHROM, POS, y = 0.5 - theta_glimpse, state, panel = panels[2])],
             bins[, .(CHROM, POS, y = 0.5 + theta_glimpse_bb, state, panel = panels[3])],
             bins[, .(CHROM, POS, y = 0.5 - theta_glimpse_bb, state, panel = panels[3])],
             pb[, .(CHROM, POS, y, state, panel = src)],
             bins[, .(CHROM, POS, y = 0.5 + theta_blocks, state, panel = panels[6])],
             bins[, .(CHROM, POS, y = 0.5 - theta_blocks, state, panel = panels[6])])
pts[, `:=`(x = offset[CHROM] + POS, panel = factor(panel, panels),
           state = factor(state, c("neu", "amp", "bamp", "del", "bdel", "loh")))]
chr <- data.table(CHROM = 1:22, start = offset, end = offset + hg38)[, mid := (start + end) / 2]
cols <- c(neu = "grey55", amp = "#d7301f", bamp = "#fc9272", del = "#2b6cb0", bdel = "#9ecae1", loh = "#2f8f46")

p <- ggplot(pts, aes(x, y)) +
  geom_rect(data = chr[CHROM %% 2 == 0], aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey95") +
  geom_hline(yintercept = c(1 / 3, 2 / 3), linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = state), size = 1.1, alpha = 0.85, stroke = 0) +
  facet_grid(panel ~ ., switch = "y") +
  scale_x_continuous(breaks = chr$mid, labels = chr$CHROM, expand = c(0, 0), position = "top") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = NULL, y = NULL, title = glue("{opt$sample} - BAF bins: 1kGP (HWE) vs GLIMPSE2 genotype prior; Eagle vs GLIMPSE phase"),
       subtitle = glue("unphased: ~{opt$bin_snps} expected hets per bin, all sites with reads, drawn at 0.5 +/- theta ",
                       "(HWE het-rate h = {h_hat}); phased: {opt$bin_snps} hets per bin. Dotted lines: 1/3, 2/3")) +
  theme_classic(base_size = 11) +
  theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing = unit(0.4, "lines"),
        axis.ticks.x = element_blank(), panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.3))
out <- file.path(sd, glue("{opt$sample}_baf_glimpse.png"))
ggsave(out, p, width = 16, height = 11.5, dpi = 150)
message("wrote ", out)
