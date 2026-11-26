#!/usr/bin/env Rscript
# Allele + expression copy-number calls per segment, arm and chromosome, from a finished numbat run.
#
# 1. Allelic imbalance: genotype-free BAF bins (R/baf_lib.R), each with a log-likelihood curve over theta.
#    RNA allele counts are overdispersed (allele-specific expression), so curves are tempered by kappa, the median
#    per-bin likelihood-ratio statistic against its chromosome's theta divided by the chi-square(1) median.
# 2. Segmentation per chromosome, jointly on theta (sum of tempered curves) and gene expression (numbat gene logFC,
#    Gaussian with a robust gene-level sd): optimal partitioning, penalty --seg_penalty per breakpoint, segments of
#    >= --min_bins bins. Genes are assigned to the nearest bin.
# 3. Copy-number state per segment, (total, minor) in {(1,0) loss, (2,0) CN-LOH, (2,1) neutral, (3,1) gain,
#    (3,0) gain+LOH, (4,2) balanced double gain, (4,1) double gain}, with blast fraction p and the expression
#    baseline fitted jointly over all segments:
#       theta = |(p (c - b) + (1 - p)) / (p c + 2 (1 - p)) - 0.5|      log2 ratio = log2((p c + 2 (1 - p)) / 2) + baseline
#    Balanced regions in RNA do not reach theta = 0 (allele-specific expression), so the observed theta is modelled as
#    sqrt(theta^2 + theta0^2), theta0 = median theta of bins below 0.12 (the balanced floor of this sample).
#    score per segment = allele term (the segment's own theta log-likelihood curve read at the predicted theta, -2 log-lik,
#    sd floored at --theta_sd) + z_expr^2 (allowing --expr_sd systematic error) + both-allele term (observed vs expected
#    rate of sites showing both alleles) + --state_penalty per copy changed (PURPLE event distance, as HAPLOID; the
#    balanced 4:2 gain counts as one with --ball_prior). p runs 0.05-1 (step 0.01) and is fitted on segments
#    >= --fit_min_mb; a segment bootstrap (--boot) gives its 95% interval, and each called segment's own implied p is
#    reported (HAPLOID clonal_tf.py) - clonal events line up on one p.
#    Expression cannot fix absolute ploidy: the best solution with a different mean copy number is reported, and
#    flagged when its score is within --ambiguity of the best (e.g. near-haploid vs doubled hypodiploid).
#    Short non-neutral segments (< --min_call_bins bins or < --min_call_mb) are labelled "focal imbalance": in RNA these
#    are mostly mono-allelic expression (imprinted clusters such as 11p15.5, protocadherins at 5q31) and RNA cannot
#    resolve focal copy number anyway; they are listed for review, not counted. Balanced double gains rest on
#    expression alone (theta unchanged), the evidence most affected by aligner / counting differences against the
#    reference (e.g. 1q21 pericentromeric genes with DRAGEN counts vs the STAR-based RNAseqCNV reference), so they need
#    2 x --min_call_bins bins; shorter ones are listed as "focal expression change".
# 4. Arm and chromosome calls (the state covering most of the arm; whole-chromosome when every arm agrees over
#    >= 80%), chromosome number, ploidy group, LOH regions and flags.
#
# Output (in <sample_dir>): calls_segments.tsv, calls_arms.tsv, calls_summary.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(glue)
})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "baf_lib.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character"),
  make_option("--sample_dir", type = "character"),
  make_option("--bin_snps", type = "integer", default = 20),
  make_option("--min_bins", type = "integer", default = 3),
  make_option("--seg_penalty", type = "double", default = 20),
  make_option("--state_penalty", type = "double", default = 3),
  make_option("--expr_sd", type = "double", default = 0.15, help = "systematic sd of segment log2 expression"),
  make_option("--ambiguity", type = "double", default = 10),
  make_option("--min_call_bins", type = "integer", default = 5),
  make_option("--min_call_mb", type = "double", default = 3),
  make_option("--doubled_flag", type = "integer", default = 5),
  make_option("--penalty", type = "character", default = "events",
              help = "state prior: events (per copy changed, PURPLE / HAPLOID) or flat (one charge per non-neutral state)"),
  make_option("--ball_prior", type = "logical", default = TRUE, help = "charge the balanced 4:2 gain as one event"),
  make_option("--p_min", type = "double", default = 0.05, help = "lowest blast fraction on the grid"),
  make_option("--fit_min_mb", type = "double", default = 5, help = "shortest segment used to fit blast fraction / baseline"),
  make_option("--theta_sd", type = "double", default = 0.02, help = "model-error floor on the allele curve's sd"),
  make_option("--both_sd", type = "double", default = 0.03, help = "systematic sd of the both-allele ratio"),
  make_option("--both_weight", type = "double", default = 1, help = "weight of the both-allele term (0 = off)"),
  make_option("--allele_source", type = "character", default = "pileup",
              help = "pileup (every panel site with reads) or vcf (HaplotypeCaller variant sites only; flagged)"),
  make_option("--boot", type = "integer", default = 200, help = "segment bootstrap resamples for the blast-fraction interval (0 = off)")
)))
sd <- opt$sample_dir

# ---- allele bins and expression per bin ----
ab <- allele_bins(sd, opt$bin_snps)
b <- ab$bins; cv <- ab$curves
b[, cfit := {
  # per-bin LR statistic against the chromosome-level theta, for the overdispersion factor
  i <- .I; tc <- which.max(colSums(cv[i, , drop = FALSE]))
  2 * (apply(cv[i, , drop = FALSE], 1, max) - cv[i, tc])
}, by = CHROM]
kappa <- max(1, median(b$cfit) / qchisq(0.5, 1))
cv <- cv / kappa
b[, cfit := NULL]

bulk_file <- file.path(sd, "numbat", "bulk_exp_only.tsv.gz")          # numbat expression model (no phasing)
if (!file.exists(bulk_file)) bulk_file <- file.path(sd, "numbat", "bulk_joint.tsv.gz")   # runs made before Nov 2026
bulk <- fread(bulk_file)
genes <- unique(bulk[!is.na(logFC) & !is.na(gene), .(CHROM, POS = (gene_start + gene_end) / 2, logFC, gene)], by = "gene")
genes <- genes[CHROM %in% 1:22]
genes[, unit := {
  bb <- b[CHROM == .BY$CHROM]
  if (nrow(bb) == 0) NA_integer_ else {
    edges <- c(-Inf, (bb$end[-nrow(bb)] + bb$start[-1]) / 2, Inf)
    findInterval(POS, edges)
  }
}, by = CHROM]
genes <- genes[!is.na(unit)]
sig <- genes[, .(r = logFC - median(logFC)), by = CHROM][, mad(r)]

# ---- segmentation ----
segment_chrom <- function(cm, x, n_g, s1, s2) {
  # cm: bins x theta curves; x: per-bin gene count, sum, sum of squares
  n <- nrow(cm); L <- min(opt$min_bins, n)
  C <- rbind(0, apply(cm, 2, cumsum)); G <- rbind(0, apply(cbind(n_g, s1, s2), 2, cumsum))
  cost <- function(i, j) {
    g <- G[j + 1, ] - G[i, ]
    sse <- if (g[1] > 0) g[3] - g[2]^2 / g[1] else 0
    -max(C[j + 1, ] - C[i, ]) + sse / (2 * sig^2)
  }
  Fv <- c(0, rep(Inf, n)); back <- integer(n)
  for (j in seq_len(n)) for (i in seq_len(max(0, j - L + 1))) {
    if (i > 1 && i - 1 < L) next                  # leaves a first segment shorter than L
    v <- Fv[i] + cost(i, j) + opt$seg_penalty
    if (v < Fv[j + 1]) { Fv[j + 1] <- v; back[j] <- i }
  }
  if (!is.finite(Fv[n + 1])) return(rep(1L, n))
  seg <- integer(n); j <- n; k <- 0L
  while (j > 0) { i <- back[j]; k <- k + 1L; seg[i:j] <- k; j <- i - 1 }
  max(seg) + 1L - seg
}
gstat <- genes[, .(n_g = .N, s1 = sum(logFC), s2 = sum(logFC^2)), by = .(CHROM, unit)]
b[, unit := seq_len(.N), by = CHROM]
b <- merge(b, gstat, by = c("CHROM", "unit"), all.x = TRUE, sort = FALSE)
setorder(b, CHROM, bin); cv <- cv[paste(b$CHROM, b$bin), , drop = FALSE]
for (v in c("n_g", "s1", "s2")) set(b, which(is.na(b[[v]])), v, 0)
b[, seg := segment_chrom(cv[.I, , drop = FALSE], NULL, n_g, s1, s2), by = CHROM]

segs <- b[, {
  f <- theta_fit(colSums(cv[.I, , drop = FALSE]))
  ng <- sum(n_g)
  .(start = min(start), end = max(end), bins = .N, n_genes = ng, theta = f$theta, theta_lo = f$lo, theta_hi = f$hi,
    theta_se = f$se, log2r = if (ng > 0) sum(s1) / ng else NA_real_, log2r_se = if (ng > 0) sig / sqrt(ng) else NA_real_)
}, by = .(CHROM, seg)]

# ---- copy-number states, blast fraction, baseline ----
states <- data.table(c = c(1, 2, 2, 3, 3, 4, 4), b = c(0, 0, 1, 1, 0, 2, 1),
                     label = c("loss", "CN-LOH", "neutral", "gain", "gain+LOH", "balanced double gain", "double gain"))
low <- b$theta[b$theta < 0.12]
theta0 <- if (length(low) >= 30) median(low) else 0.05
pred <- function(p) states[, .(k = .I, th = sqrt(((p * (c - b) + (1 - p)) / (p * c + 2 * (1 - p)) - 0.5)^2 + theta0^2),
                               er = log2((p * c + 2 * (1 - p)) / 2))]
n_seg <- nrow(segs)
seg_key <- paste(segs$CHROM, segs$seg)
ex_sd <- sqrt(fifelse(is.na(segs$log2r_se), 1, segs$log2r_se)^2 + opt$expr_sd^2)
has_ex <- !is.na(segs$log2r)
x_ex <- fifelse(has_ex, segs$log2r, 0)
seg_len <- as.numeric(segs$end - segs$start)

# value of each row of M (rows x theta_grid) at each theta in th (linear interpolation) -> rows x length(th)
at_theta <- function(M, th) {
  x <- pmin(pmax(th, 0), max(theta_grid)) / 0.005
  lo <- pmin(floor(x + 1e-9), length(theta_grid) - 1); hi <- pmin(lo + 1, length(theta_grid) - 1); fr <- x - lo
  M[, lo + 1, drop = FALSE] * rep(1 - fr, each = nrow(M)) + M[, hi + 1, drop = FALSE] * rep(fr, each = nrow(M))
}

# allele term: the segment's own log-likelihood curve over theta (sum of its tempered bin curves), read at the predicted
# theta, in -2 log-lik units. Flattened so its effective sd is at least --theta_sd (model error, as the former
# sqrt(se^2 + 0.02^2)); the curve keeps its real shape (skewed near 0 and at the 0.49 edge, flat for few hets).
Cs <- rowsum(cv, paste(b$CHROM, b$seg))[seg_key, , drop = FALSE]
c_max <- apply(Cs, 1, max)
curve_se <- apply(Cs, 1, function(cu) {
  i <- which.max(cu); d <- c(if (i > 1) cu[i] - cu[i - 1], if (i < length(cu)) cu[i] - cu[i + 1])
  0.005 / sqrt(2 * max(mean(d), 1e-6))
})
w_th <- curve_se^2 / (curve_se^2 + opt$theta_sd^2)
theta_term <- function(th) 2 * (c_max - at_theta(Cs, th)) * w_th

# both-allele term: a panel site "shows both alleles" when each allele has >= max(2, 10% of) its reads (the 0.1-0.9
# genotyping band). Expected count at minor-allele fraction 0.5 - theta, over a segment's sites (>= 4 reads):
#   sum 2hf(1-f) P(both | Bin(d, 0.5 - theta)) + (1 - 2hf(1-f)) P(both | Bin(d, eps))
# LOH turns hets into apparent homs, so the rate drops where the theta curve is least sure (one-allele sites could be
# homs). observed / expected is normalised to the balanced bins (absorbs allele-specific expression and the het rate);
# sd = phi x binomial sd (phi from balanced bins) with a floor --both_sd. Predicted ratio for a state =
# E(theta_state) / E(theta0). It reuses the same reads as the theta curve, so it is a second view, not new data.
pboth <- function(d, k, p) pbinom(d - k, d, p) - pbinom(k - 1, d, p)
st <- ab$sites[DP >= 4, .(CHROM, POS, AF, DP, ALT_N, bin)]
st[, k := pmax(2, ceiling(0.1 * DP))]
st[, both := pmin(ALT_N, DP - ALT_N) >= k]
st[, phet := pmin(2 * ab$h * AF * (1 - AF), 1)]
st[, ehom := (1 - phet) * pboth(DP, k, baf_eps)]
dps <- sort(unique(st$DP)); kk <- pmax(2, ceiling(0.1 * dps))
st[, di := match(DP, dps)]
Pd <- vapply(theta_grid, function(t) pboth(dps, kk, 0.5 - t), numeric(length(dps)))       # depth x theta
st[, p0 := phet * pboth(DP, k, 0.5 - theta0) + ehom]
st[, si := {
  sc <- segs[CHROM == .BY$CHROM]
  match(paste(.BY$CHROM, sc$seg[pmax(1L, findInterval(POS, sc$start))]), seg_key)
}, by = CHROM]
# normalise on whole balanced-looking segments (theta < 0.08, >= --min_call_bins bins), not on bins selected for low
# theta: bins picked for low theta also over-represent sites showing both alleles, which biased the reference upward
bal_seg <- which(segs$theta < 0.08 & segs$bins >= opt$min_call_bins)
if (length(bal_seg) < 3) bal_seg <- which(segs$theta < 0.12)
bal_bins <- b[paste(CHROM, seg) %in% seg_key[bal_seg], .(CHROM, bin)]
bs_tab <- st[si %in% bal_seg]
both_norm <- sum(bs_tab$both) / sum(bs_tab$p0)
bin_both <- st[!is.na(bin), .(obs = sum(both), e0 = sum(p0), v0 = sum(p0 * (1 - p0))), by = .(CHROM, bin)][
  bal_bins, on = .(CHROM, bin), nomatch = NULL]
both_phi <- max(1, mad((bin_both$obs - both_norm * bin_both$e0) / sqrt(bin_both$v0)))
Eseg <- matrix(NA_real_, n_seg, length(theta_grid))
agg <- st[!is.na(si), .(w = sum(phet)), by = .(si, di)]
Ehet <- rowsum(agg$w * Pd[agg$di, , drop = FALSE], agg$si)
eh <- st[!is.na(si), .(eh = sum(ehom), obs = sum(both), e0 = sum(p0), v0 = sum(p0 * (1 - p0))), keyby = si]
Eseg[as.integer(rownames(Ehet)), ] <- Ehet
Eseg[eh$si, ] <- Eseg[eh$si, , drop = FALSE] + eh$eh
both_obs <- both_se <- rep(NA_real_, n_seg)
both_obs[eh$si] <- eh$obs / (both_norm * eh$e0)
both_se[eh$si] <- sqrt(both_phi^2 * eh$v0 / (both_norm * eh$e0)^2 + opt$both_sd^2)
has_both <- !is.na(both_obs) & opt$both_weight > 0
E_ref <- as.vector(at_theta(Eseg, theta0))
both_pred <- function(th) at_theta(Eseg, th) / E_ref
both_term <- function(th) {
  z <- opt$both_weight * (both_obs - both_pred(th))^2 / both_se^2
  z[!has_both, ] <- 0; z
}

# grid over blast fraction p and expression baseline. Only segments of >= --min_call_bins bins and >= --fit_min_mb
# (5 Mb, as HAPLOID's MIN_SEG) are fitted: shorter ones are mostly allele-specific expression, and one of them can pull
# p on its own (on HH-1 a 3.1 Mb segment at 11q12-13, theta 0.145, fits a loss at p ~ 0.48 and nothing at p = 1).
# Every segment then gets its best state at the fitted (p, baseline).
# state prior: --penalty events (default) charges --state_penalty per copy changed, the PURPLE event distance HAPLOID
# uses (|major - 1| + |minor - 1|: gain 3:1 and loss 1:0 = 1, CN-LOH 2:0, 4:2 and 4:1 = 2, 3:0 = 3). Without it a
# one-copy gain at p = 1 and a two-copy gain (4:1) at p ~ 0.5 have the same theta and expression and the same cost, so
# p is held up by whichever LOH or 2:2 segment happens to be present. --penalty flat = the former one charge per state.
# B-ALL prior (as --ball_prior in run_numbat.R): the balanced 4:2 gain (tetrasomy 21, the hallmark of high
# hyperdiploidy) costs one event like a single-copy gain; at 2 its expression-only evidence lost +21x2 below p ~ 0.4.
ev <- with(states, abs(pmax(c - b, b) - 1) + abs(pmin(c - b, b) - 1))
if (opt$ball_prior) ev[states$label == "balanced double gain"] <- 1
pen <- matrix(opt$state_penalty * (if (opt$penalty == "events") ev else as.numeric(ev > 0)), n_seg, nrow(states),
              byrow = TRUE)
fit_rows <- which(segs$bins >= opt$min_call_bins & seg_len >= opt$fit_min_mb * 1e6)
if (length(fit_rows) < 3) fit_rows <- seq_len(n_seg)
p_grid <- seq(opt$p_min, 1, by = 0.01); base_grid <- seq(-1.5, 1, by = 0.025)
z_at <- function(A, er, base) A + has_ex * outer(x_ex - base, er, "-")^2 / ex_sd^2      # has_ex recycles down rows
grid <- CJ(pi = seq_along(p_grid), bi = seq_along(base_grid))
grid[, p := p_grid[pi]][, base := base_grid[bi]]
Zseg <- matrix(NA_real_, length(fit_rows), nrow(grid))   # per-segment best score at every grid cell (bootstrap)
mcn <- numeric(nrow(grid))
L_fit <- seg_len[fit_rows]
cell <- 0L
for (i in seq_along(p_grid)) {
  pr <- pred(p_grid[i])
  A <- (theta_term(pr$th) + both_term(pr$th) + pen)[fit_rows, , drop = FALSE]
  exf <- has_ex[fit_rows]; xf <- x_ex[fit_rows]; sdf <- ex_sd[fit_rows]
  for (j in seq_along(base_grid)) {
    cell <- cell + 1L
    z <- A + exf * outer(xf - base_grid[j], pr$er, "-")^2 / sdf^2
    k <- max.col(-z, ties.method = "first")
    Zseg[, cell] <- z[cbind(seq_along(k), k)]
    mcn[cell] <- sum(states$c[k] * L_fit) / sum(L_fit)
  }
}
grid[, `:=`(score = colSums(Zseg), mean_cn = mcn)]
best <- which.min(grid$score)
p_hat <- grid$p[best]; base_hat <- grid$base[best]
alt <- grid[abs(mean_cn - grid$mean_cn[best]) >= 0.5][which.min(score)]
pr <- pred(p_hat)
z_all <- z_at(theta_term(pr$th) + both_term(pr$th) + pen, pr$er, base_hat)
fit <- list(k = max.col(-z_all, ties.method = "first"), mean_cn = grid$mean_cn[best])
segs[, `:=`(c = states$c[fit$k], minor = states$b[fit$k], call = states$label[fit$k])]
segs[, `:=`(theta_exp = round(pr$th[fit$k], 3), log2r_exp = round(pr$er[fit$k] + base_hat, 3),
            both_ratio = round(both_obs, 3), both_ratio_exp = round(both_pred(pr$th)[cbind(seq_len(n_seg), fit$k)], 3))]
segs[call == "balanced double gain" & bins < 2 * opt$min_call_bins, `:=`(call = "focal expression change", c = NA, minor = NA)]
segs[!call %in% c("neutral", "focal expression change") & (bins < opt$min_call_bins | end - start < opt$min_call_mb * 1e6),
     `:=`(call = "focal imbalance", c = NA, minor = NA)]
imbalanced <- any(!segs$call %in% c("neutral", "balanced double gain", "focal imbalance", "focal expression change"))

# segment bootstrap of p (HAPLOID joint_fit.py): resample the fitted segments with replacement and refind the best grid
# cell on the summed per-segment surfaces. The interval shows how much p depends on which events drive it.
p_lo <- p_hi <- NA_real_
if (opt$boot > 0 && length(fit_rows) >= 3) {
  set.seed(12345)
  W <- rmultinom(opt$boot, length(fit_rows), rep(1, length(fit_rows)))       # fitted segments x resamples
  kb <- max.col(-(t(W) %*% Zseg), ties.method = "first")
  pb <- grid$p[kb]
  p_lo <- round(unname(quantile(pb, 0.025)), 2); p_hi <- round(unname(quantile(pb, 0.975)), 2)
  fwrite(data.table(resample = seq_along(pb), p = pb, mean_cn = round(grid$mean_cn[kb], 2),
                    segments_left_out = apply(W, 2, function(w) paste(seg_key[fit_rows][w == 0], collapse = ","))),
         file.path(sd, "blast_fraction_boot.tsv"), sep = "\t")
}
profile <- grid[, .(score = min(score)), by = p][, d_score := score - min(score)]
if (!imbalanced) p_hat <- p_lo <- p_hi <- NA_real_

# blast fraction implied by each called segment on its own (HAPLOID clonal_tf.py): invert the expected theta of its
# state (expression for balanced 2:2 gains, whose theta does not move with p), with the 1-SE support interval. Clonal
# events line up on one p; outliers are subclonal, germline (constitutional / homozygous runs) or misassigned.
pf <- seq(0, 1, by = 0.001)
implied <- function(cc, mm, th, lo, hi, lr, lr_se) {
  if (is.na(cc) || (cc == 2 && mm == 1)) return(list(NA_real_, NA_real_, NA_real_, NA_character_))
  if (cc - 2 * mm == 0) {                              # balanced (4:2): expression only
    er <- log2((pf * cc + 2 * (1 - pf)) / 2); y <- lr - base_hat
    inv <- function(v) approx(er, pf, xout = v, rule = 2, ties = "ordered")$y
    return(list(inv(y), inv(min(y - lr_se, y + lr_se)), inv(max(y - lr_se, y + lr_se)), "expression"))
  }
  tf <- sqrt(((pf * (cc - mm) + (1 - pf)) / (pf * cc + 2 * (1 - pf)) - 0.5)^2 + theta0^2)
  inv <- function(v) approx(tf, pf, xout = v, rule = 2, ties = "ordered")$y
  list(inv(th), inv(lo), inv(hi), "BAF")
}
imp <- rbindlist(lapply(seq_len(n_seg), function(i) {
  r <- implied(segs$c[i], segs$minor[i], segs$theta[i], segs$theta_lo[i], segs$theta_hi[i], segs$log2r[i], segs$log2r_se[i])
  data.table(p_implied = r[[1]], p_implied_lo = r[[2]], p_implied_hi = r[[3]], p_source = r[[4]])
}))
segs[, names(imp) := lapply(imp, function(v) if (is.numeric(v)) round(v, 3) else v)]

# ---- arms and chromosomes ----
arm_rows <- rbindlist(lapply(1:22, function(ch) {
  cen <- hg38_cen[ch]
  arms <- if (ch %in% acrocentric) list(q = c(cen, hg38_len[ch])) else list(p = c(0, cen), q = c(cen, hg38_len[ch]))
  rbindlist(lapply(names(arms), function(a) {
    r <- arms[[a]]; bi <- which(b$CHROM == ch & b$POS >= r[1] & b$POS < r[2])
    if (!length(bi)) return(NULL)
    s <- segs[CHROM == ch & end >= r[1] & start < r[2] & !call %like% "^focal"]   # focal calls are not counted
    s[, ov := pmin(end, r[2]) - pmax(start, r[1])]
    cov <- s[, .(ov = sum(ov)), by = call][order(-ov)]
    f <- theta_fit(colSums(cv[bi, , drop = FALSE]))
    g <- genes[CHROM == ch & POS >= r[1] & POS < r[2]]
    data.table(CHROM = ch, arm = a, start = min(b$start[bi]), end = max(b$end[bi]), bins = length(bi),
               n_genes = nrow(g), theta = f$theta, theta_se = round(f$se, 3),
               log2r = round(mean(g$logFC) - base_hat, 3), call = cov$call[1],
               call_frac = round(cov$ov[1] / sum(cov$ov), 2),
               c = s[call == cov$call[1], c[1]], minor = s[call == cov$call[1], minor[1]],
               events = paste(s[call != "neutral", sprintf("%s:%.1f-%.1fMb", call, start / 1e6, end / 1e6)],
                              collapse = ";"))
  }))
}))
chrom <- arm_rows[, .(whole = all(call == call[1]) && all(call_frac >= 0.8), call = call[1], c = c[1], minor = minor[1]),
                  by = CHROM]
nb <- tryCatch(fread(file.path(sd, "karyotype.tsv"))[, .(CHROM, numbat_call = call)], error = function(e) NULL)
if (!is.null(nb)) arm_rows <- merge(arm_rows, nb, by = "CHROM", all.x = TRUE)

# ---- summary ----
wc <- chrom[whole == TRUE]
n_chrom <- 46 + wc[, sum(c - 2)]
group <- as.character(cut(n_chrom, c(-Inf, 23, 30, 39, 43, 45, 46, 50, 67, 79, Inf),
  c("<24", "near-haploid (24-30)", "low hypodiploid (31-39)", "high hypodiploid (40-43)",
    "near-diploid hypodiploid (44-45)", "diploid (46)", "low hyperdiploid (47-50)", "high hyperdiploid (51-67)",
    "near-triploid (68-79)", ">79")))
fmt <- function(ch, times, sign) { x <- rep(ch, times); if (length(x)) paste0(sign, sort(x), collapse = ",") else "" }
loh <- segs[minor %in% 0]
loh_regions <- loh[!(CHROM %in% wc[minor == 0, CHROM]),
                   paste(sprintf("%d:%.1f-%.1fMb(%s)", CHROM, start / 1e6, end / 1e6, call), collapse = ";")]
n_whole_loh <- wc[minor == 0 & c >= 2, .N]
ambig <- nrow(alt) == 1 && alt$score - grid$score[best] < opt$ambiguity
partial <- chrom[whole == FALSE, CHROM]
numbat_n <- tryCatch(fread(file.path(sd, "summary.tsv"))$est_chromosomes, error = function(e) NA)
flags <- c(
  if (n_whole_loh >= opt$doubled_flag) glue("possible_doubled_hypodiploid({n_whole_loh} whole-chromosome LOH)"),
  if (ambig) glue("ploidy_ambiguous(alternative mean copy number {round(alt$mean_cn, 2)} vs {round(fit$mean_cn, 2)}, ",
                  "score +{round(alt$score - grid$score[best], 1)})"),
  if (nrow(loh) > 0) "LOH_regions_present(somatic or germline homozygosity: confirm with remission DNA)",
  if (length(partial)) glue("partial_chromosome_events({paste(partial, collapse = ',')})"),
  if (!is.na(numbat_n) && numbat_n != n_chrom) glue("differs_from_numbat({numbat_n})"),
  if (opt$allele_source == "vcf") "vcf_input(variant sites only: LOH where the reference allele is kept can be missed or shortened; het rate h is not meaningful)"
)
summary <- data.table(
  sample = opt$sample, est_chromosomes = n_chrom, ploidy_group = group,
  gains = fmt(wc[c > 2, CHROM], wc[c > 2, c - 2], "+"), losses = fmt(wc[c < 2, CHROM], wc[c < 2, 2 - c], "-"),
  whole_chrom_loh = fmt(wc[minor == 0, CHROM], 1, "LOH"), loh_regions = loh_regions,
  focal = segs[call %like% "^focal", paste(sprintf("%d:%.1f-%.1fMb(%s)", CHROM, start / 1e6, end / 1e6,
                                                  sub("focal ", "", call)), collapse = ";")],
  blast_fraction = p_hat, blast_fraction_lo = p_lo, blast_fraction_hi = p_hi, fit_segments = length(fit_rows),
  het_rate_h = ab$h, overdispersion_kappa = round(kappa, 2), theta0 = theta0, both_allele_phi = round(both_phi, 2),
  expr_baseline = base_hat,
  alt_mean_cn = if (nrow(alt)) round(alt$mean_cn, 2) else NA, alt_score_gap = if (nrow(alt)) round(alt$score - grid$score[best], 1) else NA,
  flags = paste(flags, collapse = ";"),
  note = "autosomes only; blast fraction from imbalanced segments; CN-LOH may be germline homozygosity"
)
fwrite(segs[, .(CHROM, start, end, bins, n_genes, theta, theta_lo, theta_hi, theta_exp, log2r = round(log2r - base_hat, 3),
                log2r_exp = round(log2r_exp - base_hat, 3), c, minor, call, both_ratio, both_ratio_exp,
                p_implied, p_implied_lo, p_implied_hi, p_source)], file.path(sd, "calls_segments.tsv"), sep = "\t")
fwrite(profile[, .(p, d_score = round(d_score, 2))], file.path(sd, "blast_fraction_profile.tsv"), sep = "\t")

# ---- blast-fraction plot: per-segment implied p (left) and the score profile over p (right) ----
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
ps <- segs[!is.na(p_implied)][, lab := sprintf("%d:%.0f-%.0f Mb", CHROM, start / 1e6, end / 1e6)]
ps[, lab := factor(lab, levels = lab[order(CHROM, start)])]
band <- if (!is.na(p_lo)) annotate("rect", xmin = -Inf, xmax = Inf, ymin = p_lo, ymax = p_hi, alpha = 0.15, fill = "#3060a0")
g1 <- ggplot(ps, aes(lab, p_implied, colour = call, shape = p_source)) + band +
  { if (!is.na(p_hat)) geom_hline(yintercept = p_hat, colour = "#3060a0", linewidth = 0.6) } +
  geom_errorbar(aes(ymin = p_implied_lo, ymax = p_implied_hi), width = 0.2) + geom_point(size = 2.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c(gain = "#c0392b", `CN-LOH` = "#2e8b57", loss = "#2c6fbb", `gain+LOH` = "#8e44ad",
                                 `balanced double gain` = "#e67e5a", `double gain` = "#7b241c"), drop = TRUE) +
  labs(x = NULL, y = "implied blast fraction", colour = NULL, shape = NULL,
       title = "Blast fraction implied by each called segment",
       subtitle = "bars: 1-SE support; line and band: joint fit and 95% segment bootstrap") +
  theme_bw(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
g2 <- ggplot(profile, aes(p, d_score)) + geom_line(colour = "grey25") +
  { if (!is.na(p_lo)) annotate("rect", xmin = p_lo, xmax = p_hi, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#3060a0") } +
  coord_cartesian(ylim = c(0, min(200, max(profile$d_score)))) +
  labs(x = "blast fraction p", y = "score - best (baseline profiled out)", title = "Fit score over p") + theme_bw(11)
ttl <- if (is.na(p_hat)) "blast fraction: not estimable (no imbalanced segment)" else
  sprintf("blast fraction %.2f (95%% bootstrap %.2f-%.2f; %d segments fitted)", p_hat, p_lo, p_hi, length(fit_rows))
ggsave(file.path(sd, glue("{opt$sample}_blast_fraction.png")),
       (g1 | g2) + plot_layout(widths = c(1.6, 1)) + plot_annotation(title = paste(opt$sample, "-", ttl)),
       width = 13, height = 5.5, dpi = 150)
fwrite(arm_rows, file.path(sd, "calls_arms.tsv"), sep = "\t")
fwrite(summary, file.path(sd, "calls_summary.tsv"), sep = "\t")
message(sprintf("h = %.2f, kappa = %.2f, theta0 = %.3f, %d segments, blast fraction %s, baseline %.3f",
                ab$h, kappa, theta0, nrow(segs), format(p_hat), base_hat))
print(summary[, .(sample, est_chromosomes, ploidy_group, gains, losses, whole_chrom_loh, loh_regions, blast_fraction, flags)])
