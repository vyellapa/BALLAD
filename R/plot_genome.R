#!/usr/bin/env Rscript
# Genome-wide copy-number plot from a finished sample directory, coloured by the copy-number call of the segment each
# point falls in (calls_segments.tsv from R/call_arms.R; numbat's joint segments when that is absent).
#
#   Expression  gene-level log2 expression ratio vs the B-ALL reference (minus the fitted baseline when calls exist).
#               Faint dots: single genes. Coloured dots: medians of --bin_genes consecutive genes. Lines: segment means.
#   BAF         Faint dots: ALT / (REF + ALT) at single panel sites with depth >= --min_dp (homozygous sites sit at
#               0 / 1, hets around 0.5). Coloured dots: genotype-free allelic-imbalance bins (R/baf_lib.R: every site with
#               reads, any depth, 1kGP prior, no phasing), drawn at 0.5 +/- theta; lines = segment theta. A clonal
#               single-copy gain sits at 1/3 and 2/3, a balanced gain (2:2) at 0.5, LOH near 0 and 1. Balanced regions
#               sit slightly off 0.5 (RNA allele-specific expression floor, theta0 in calls_summary.tsv).
# Called events are labelled above the plot with their copy change (+1, +2, -1, CN-LOH), so calls never rest on colour
# alone. Colours (validated): deep red gain, ochre two-copy balanced gain, navy loss, deep teal CN-LOH, cool grey neutral. Focal segments (not counted) are left out (bins and lines). The
# subtitle has the overall tumor fraction and its 95% interval.
#
# Output: <sample_dir>/<sample>_genome.png

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(glue)
})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "baf_lib.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character"),
  make_option("--sample_dir", type = "character"),
  make_option("--min_dp", type = "integer", default = 10),
  make_option("--exp_limit", type = "double", default = 1.5),
  make_option("--bin_snps", type = "integer", default = 20),
  make_option("--bin_genes", type = "integer", default = 25),
  make_option("--out", type = "character", default = NULL, help = "default <sample_dir>/<sample>_genome.png")
)))
sd <- opt$sample_dir
offset <- c(0, cumsum(as.numeric(hg38_len)))[1:22]
gpos <- function(chrom, pos) offset[chrom] + pos

# ---- segments and states: our calls when present, else numbat ----
calls_file <- file.path(sd, "calls_segments.tsv")
has_calls <- file.exists(calls_file)
has <- function(x) length(x) && !is.na(x) && nzchar(as.character(x))
if (has_calls) {
  cs <- fread(calls_file)
  sm <- fread(file.path(sd, "calls_summary.tsv"))
  base <- sm$expr_baseline
  cs[, state := fcase(call == "neutral", "neu", call == "balanced double gain", "bamp",
                      call %in% c("gain", "double gain", "gain+LOH"), "amp", call == "loss", "del",
                      call == "CN-LOH", "loh", default = "focal")]
  segs <- cs[, .(CHROM, start, end, state, theta, log2r, cn = c)]
  title_calls <- paste0(sm$est_chromosomes, " chromosomes  \u00b7  ", sub(" \\(.*\\)$", "", sm$ploidy_group))
  bf <- if ("blast_fraction" %in% names(sm) && has(sm$blast_fraction)) {
    if ("blast_fraction_lo" %in% names(sm) && has(sm$blast_fraction_lo))
      sprintf("tumor fraction %.2f, 95%% CI %.2f-%.2f", sm$blast_fraction, sm$blast_fraction_lo, sm$blast_fraction_hi)
    else sprintf("tumor fraction %.2f", sm$blast_fraction)
  } else NULL
  title_calls <- paste(c(title_calls, bf), collapse = "  ·  ")
} else {
  nf <- file.path(sd, "numbat", "segments_exp_only.tsv")
  if (!file.exists(nf)) nf <- file.path(sd, "numbat", "segments_joint.tsv")
  ns <- fread(nf)
  base <- 0
  segs <- ns[, .(CHROM, start = seg_start, end = seg_end,
                 state = fcase(sub("_.*", "", cnv_state_post) == "amp", "amp", cnv_state_post == "bamp", "bamp",
                               sub("_.*", "", cnv_state_post) %in% c("del", "bdel"), "del", sub("_.*", "", cnv_state_post) == "loh", "loh", default = "neu"),
                 theta = NA_real_, log2r = log2(phi_mle))]
  segs[, cn := fcase(state == "amp", 3, state == "bamp", 4, state == "del", 1, default = 2)]
  title_calls <- "numbat calls"
}
state_of <- function(d) {
  s <- segs[, .(CHROM, start, end, state)]; setkey(s, CHROM, start, end)
  d <- copy(d)[, `:=`(start = POS, end = POS)]
  o <- foverlaps(d, s, by.x = c("CHROM", "start", "end"), mult = "first", nomatch = NA)
  o[, state := fifelse(is.na(state), "neu", state)][, .(CHROM, POS, y, state)]
}
clamp <- function(x) pmax(pmin(x, opt$exp_limit), -opt$exp_limit)

# ---- data ----
bulk_file <- file.path(sd, "numbat", "bulk_exp_only.tsv.gz")          # numbat expression model (no phasing)
if (!file.exists(bulk_file)) bulk_file <- file.path(sd, "numbat", "bulk_joint.tsv.gz")   # runs made before Nov 2026
bulk <- fread(bulk_file)
genes <- unique(bulk[!is.na(logFC) & !is.na(gene), .(CHROM, POS = gene_start, logFC = logFC - base, gene)],
                by = "gene")[CHROM %in% 1:22][order(CHROM, POS)]
ex <- genes[, .(CHROM, POS, y = clamp(logFC))]
ex_bins <- genes[, bin := (seq_len(.N) - 1) %/% opt$bin_genes, by = CHROM][
  , if (.N >= opt$bin_genes / 2) .(POS = as.numeric(median(POS)), y = clamp(median(logFC))), by = .(CHROM, bin)]
ab <- allele_bins(sd, opt$bin_snps)
bins <- ab$bins
baf <- ab$sites[DP >= opt$min_dp, .(CHROM, POS, y = ALT_N / DP)]

panels <- c("Expression (log2 ratio)", "Allele fraction (BAF)")
lv <- c("neu", "amp", "bamp", "del", "loh", "focal")
pts <- rbind(ex[, panel := panels[1]], baf[, panel := panels[2]])
pts[, `:=`(x = gpos(CHROM, POS), panel = factor(panel, panels))]
bpts <- rbind(state_of(ex_bins[, .(CHROM, POS, y)])[, panel := panels[1]],
              state_of(bins[, .(CHROM, POS, y = 0.5 + theta)])[, panel := panels[2]],
              state_of(bins[, .(CHROM, POS, y = 0.5 - theta)])[, panel := panels[2]])
bpts[, `:=`(x = gpos(CHROM, POS), panel = factor(panel, panels), state = factor(state, lv))]
bpts <- bpts[state != "focal"]                                # focal segments are not counted: not drawn
setorder(bpts, state)                                         # neutral first, so calls draw on top
segs_l <- segs[state != "focal"]

lines <- rbind(segs_l[!is.na(log2r), .(x = gpos(CHROM, start), xend = gpos(CHROM, end), y = clamp(log2r), panel = panels[1])],
               segs_l[!is.na(theta), .(x = gpos(CHROM, start), xend = gpos(CHROM, end), y = 0.5 + theta, panel = panels[2])],
               segs_l[!is.na(theta), .(x = gpos(CHROM, start), xend = gpos(CHROM, end), y = 0.5 - theta, panel = panels[2])])
lines[, panel := factor(panel, panels)]
ref <- data.table(panel = factor(panels, panels), y = c(0, 0.5))
thirds <- data.table(panel = factor(panels[2], panels), y = c(1 / 3, 2 / 3), lab = c("1/3", "2/3"))
lims <- data.table(x = 0, y = c(-opt$exp_limit, opt$exp_limit, 0, 1), panel = factor(panels[c(1, 1, 2, 2)], panels))
chr <- data.table(CHROM = 1:22, start = offset, end = offset + hg38_len)[, mid := (start + end) / 2]

# ---- event labels (secondary encoding: every call is named in text) ----
ev <- segs_l[state != "neu"]
seg_len <- segs[, .(tot = sum(end - start)), by = CHROM]           # the part of each chromosome the data cover
ev_lab <- if (nrow(ev)) ev[, {
  cov <- sum(end - start) / seg_len[.BY, on = "CHROM", tot]          # acrocentric p arms have no genes: not counted
  arm <- if (all(end <= hg38_cen[CHROM])) "p" else if (all(start >= hg38_cen[CHROM])) "q" else ""
  # label = copy change (the chromosome is on the axis): +1 gain, +2 two-copy gain (2:2 or 4:1), -1 loss, CN-LOH;
  # events covering < 80% of the chromosome (of its segmented length) add the arm on a second line
  d <- cn[which.max(end - start)] - 2
  txt <- if (state[1] == "loh") "CN-LOH" else if (d > 0) paste0("+", d) else sub("-", "\u2212", as.character(d))
  if (cov < 0.8) txt <- paste0(txt, "\n", CHROM, arm)
  .(x = gpos(CHROM, (min(start) + max(end)) / 2), lab = txt)
}, by = .(CHROM, state)] else data.table(x = numeric(), lab = character())
if (nrow(ev_lab)) ev_lab[, panel := factor(panels[1], panels)][, y := opt$exp_limit * 0.93]

# ---- palette: gain #C0392B and loss #1D4E89 chosen by the user; ochre and deep teal picked to match. Validated with
# the dataviz skill's validate_palette.js (light surface, all pairs): worst CVD separation 11.1, normal-vision 18.7;
# the navy sits just below the lightness band (OKLCH L 0.423 vs 0.43), harmless on white ----
ink <- "#1f2430"; ink2 <- "#4a4f5c"; muted <- "#8b909c"; band <- "#f4f5f8"; raw <- "#cfd2d9"; seg_ink <- "#262b36"
cols <- c(neu = "#a3a7b1", amp = "#C0392B", bamp = "#D68910", del = "#1D4E89", loh = "#16A085")
labs <- c(neu = "neutral", amp = "gain", bamp = "balanced gain (2:2)", del = "loss", loh = "CN-LOH")
present <- intersect(lv, unique(as.character(bpts$state)))

p <- ggplot(pts, aes(x, y)) +
  geom_rect(data = chr[CHROM %% 2 == 0], aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = band) +
  geom_point(color = raw, size = 0.35, alpha = 0.35, stroke = 0) +
  geom_blank(data = lims) +
  geom_hline(data = ref, aes(yintercept = y), color = muted, linewidth = 0.3) +
  geom_hline(data = thirds, aes(yintercept = y), color = muted, linewidth = 0.25, linetype = "22") +
  geom_text(data = thirds, aes(x = max(chr$end), y = y, label = lab), inherit.aes = FALSE, hjust = -0.25, size = 2.8,
            color = ink2) +
  geom_point(data = bpts, aes(color = state), size = 1.3, alpha = 0.85, stroke = 0) +
  geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = y), color = seg_ink, linewidth = 0.55, lineend = "round") +
  { if (nrow(ev_lab)) geom_label(data = ev_lab, aes(x = x, y = y, label = lab), inherit.aes = FALSE, size = 3.3, lineheight = 0.9,
                                 fontface = "bold",
                                 color = ink, fill = "white", linewidth = 0, label.padding = unit(0.12, "lines")) } +
  facet_grid(panel ~ ., scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = chr$mid, labels = chr$CHROM, expand = c(0, 0), position = "top") +
  scale_color_manual(values = cols[present], labels = labs[present], breaks = present, name = NULL) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = opt$sample, subtitle = title_calls) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", color = ink, size = 13),
        plot.subtitle = element_text(color = ink2, margin = margin(b = 8)),
        plot.title.position = "plot",
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90, color = ink2, size = 10),
        axis.text.x.top = element_text(color = ink2, size = 9), axis.text.y = element_text(color = muted, size = 8),
        panel.grid = element_blank(), panel.spacing = unit(0.9, "lines"),
        legend.position = "bottom", legend.text = element_text(color = ink2, size = 9),
        legend.margin = margin(t = -4), plot.margin = margin(12, 24, 8, 8),
        plot.background = element_rect(fill = "white", color = NA))

out <- if (is.null(opt$out)) file.path(sd, glue("{opt$sample}_genome.png")) else opt$out
ggsave(out, p, width = 15, height = 6.8, dpi = 160)
message("wrote ", out)
