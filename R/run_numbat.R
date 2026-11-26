#!/usr/bin/env Rscript
# numbat bulk expression model for one B-ALL RNA-seq sample: gene-level logFC against the B-ALL reference (read by
# R/call_arms.R and R/plot_genome.R) and numbat's own expression-only calls, kept as a cross-check (karyotype.tsv,
# summary.tsv). The copy-number calls themselves come from R/call_arms.R.
#
# Inputs (written by ballad.py): allele_counts.tsv.gz (panel SNPs with reads: CHROM POS REF ALT AF REF_N ALT_N DP OTH),
# HTSeq gene counts, and the B-ALL reference (build_reference.R).
#
# No phasing. numbat's joint (expression + phased allele) mode needed Eagle2-phased hets; it was dropped with the phasing
# step, because RNA phase holds only within ~50 kb and the calls (call_arms.R, genotype-free allele model) never used it.
# analyze_bulk(exp_only = TRUE) treats every region as diploid baseline and ignores the allele columns, but get_bulk /
# analyze_bulk still require some allele rows, so hets by numbat's genotype rule (0.1 <= ALT/DP <= 0.9) are passed with a
# placeholder phase (1|0) and genetic position (1 cM/Mb). On HH-1 this gives segments and logFC identical to the
# former Eagle-phased exp_only run (max |logFC difference| 5e-15). Expression alone cannot see CN-LOH or tell a 2:2 gain
# from a 3:1, and it called -19 / -22 there that the allele model rejects: read numbat's calls as a check only.
#
# Karyotype (karyotype.tsv, one row per autosome): the numbat state covering most of the chromosome, copy change
# (amp +1, bamp +2, del -1, bdel -2, loh 0), and the both-allele rate (see below). summary.tsv: estimated chromosome
# number = 46 + sum of whole-chromosome changes (sex chromosomes are not assessed by numbat and are assumed 2),
# the ploidy group, and flags.
# Both-allele rate (LOH evidence, lp_loh het_evidence_scan): at panel sites with DP >= 10, observed sites with
# >= 2 reads of each allele vs expected sum 2pq * P(>= 2 of each | DP, 0.5), p = 1kGP AF. A chromosome whose hets have
# dropped out (LOH at high blast %) has a low ratio; normalised to the median of chromosomes called neutral.

suppressPackageStartupMessages({
  library(optparse)
  library(numbat)
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(glue)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character"),
  make_option("--counts", type = "character"),
  make_option("--sample_dir", type = "character", help = "pipeline output dir for this sample"),
  make_option("--reference", type = "character"),
  make_option("--t", type = "double", default = 1e-5),
  make_option("--cnloh_flag", type = "integer", default = 5L,
              help = "whole-chromosome CN-LOH count at which to flag possible doubled hypodiploidy"),
  make_option("--allele_source", type = "character", default = "pileup",
              help = "pileup | vcf; with vcf (variant sites only) the both-allele rate is not computed")
)))
sample <- opt$sample
sd <- opt$sample_dir
nb <- file.path(sd, "numbat"); dir.create(nb, showWarnings = FALSE)

sym_file <- file.path(dirname(opt$reference), "ens_to_symbol.tsv")  # written by build_reference.R
ens_to_symbol <- if (file.exists(sym_file)) {
  sym <- fread(sym_file)
  function(ensg) sym$symbol[match(sub("\\..*", "", ensg), sym$ENSG)]
} else {
  function(ensg) suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = sub("\\..*", "", ensg),
                                                        keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
}

# ---- expression ----
ht <- fread(opt$counts, header = FALSE, col.names = c("ENSG", "count"))[!startsWith(ENSG, "__")]
ht[, symbol := ens_to_symbol(ENSG)]
ht <- ht[!is.na(symbol) & symbol != "", .(count = sum(count)), by = symbol]
count_mat <- Matrix(ht$count, ncol = 1, sparse = TRUE, dimnames = list(ht$symbol, sample))
ref <- readRDS(opt$reference)

# ---- alleles ----
ac <- fread(file.path(sd, "allele_counts.tsv.gz"))
# placeholder allele table (see header): numbat's het rule, phase and genetic position unused in exp_only mode
hets <- ac[DP > 0 & ALT_N / DP >= 0.1 & ALT_N / DP <= 0.9][, CHROM := as.integer(sub("^chr", "", CHROM))][CHROM %in% 1:22]
g <- as.data.table(gtf_hg38)[, .(CHROM = as.integer(as.character(CHROM)), gene_start, gene_end, gene)]
setkey(g, CHROM, gene_start, gene_end)
hets[, `:=`(s = POS, e = POS)]
df_allele <- foverlaps(hets, g, by.x = c("CHROM", "s", "e"), mult = "first")[
  , .(cell = sample, snp_id = paste(CHROM, POS, REF, ALT, sep = "_"), CHROM = factor(CHROM), POS, cM = POS / 1e6,
      REF, ALT, AD = ALT_N, DP, GT = "1|0", gene)]
df_allele <- as.data.frame(df_allele)

# ---- numbat ----
segs_of <- function(b) {
  if (!"cnv_state_post" %in% names(b)) b <- mutate(b, cnv_state_post = cnv_state)   # exp_only skips the retest
  b <- mutate(b, cnv_state_post = sub("_[0-9]+$", "", cnv_state_post))                # exp_only states: amp_1 -> amp
  per_seg <- intersect(c("theta_mle", "phi_mle", "p_amp", "p_del", "p_loh", "p_bamp", "p_neu", "LLR"), names(b))
  b %>% filter(!is.na(seg)) %>%
    group_by(CHROM, seg, seg_start, seg_end, cnv_state_post) %>%
    summarise(n_genes = n_distinct(gene, na.rm = TRUE),
              logFC = mean(logFC, na.rm = TRUE), across(all_of(per_seg), dplyr::first), .groups = "drop") %>%
    arrange(CHROM, seg_start) %>% mutate(CHROM = as.integer(as.character(CHROM)))
}

bulk <- get_bulk(count_mat = count_mat, lambdas_ref = ref$lambdas, df_allele = df_allele, gtf = gtf_hg38, nu = 1,
                 verbose = FALSE)
res <- analyze_bulk(bulk, t = opt$t, find_diploid = TRUE, exp_only = TRUE, retest = TRUE, verbose = FALSE)
segs <- segs_of(res)
# keep gene rows only and drop the placeholder allele columns: they are not data
fwrite(res %>% filter(!is.na(Y_obs)) %>% dplyr::select(-any_of(c("snp_id", "REF", "ALT", "cM", "GT", "AD", "DP", "AR",
                                                                  "pAD", "pBAF", "inter_snp_cm", "p_s"))),
       file.path(nb, "bulk_exp_only.tsv.gz"), sep = "\t")
fwrite(segs, file.path(nb, "segments_exp_only.tsv"), sep = "\t")

# ---- karyotype ----
delta <- c(neu = 0, amp = 1, bamp = 2, del = -1, bdel = -2, loh = 0)
call_chrom <- function(segs) {
  segs %>% mutate(len = seg_end - seg_start) %>%
    group_by(CHROM) %>%
    summarise(
      call = names(which.max(tapply(len, cnv_state_post, sum))),
      call_frac = round(max(tapply(len, cnv_state_post, sum)) / sum(len), 2),
      events = paste(ifelse(cnv_state_post == "neu", NA, glue("{cnv_state_post}:{round(seg_start/1e6, 1)}-{round(seg_end/1e6, 1)}Mb")) %>%
                       na.omit(), collapse = ";"),
      phi = round(phi_mle[which.max(len)], 2),
      theta = if ("theta_mle" %in% names(segs)) round(theta_mle[which.max(len)], 3) else NA_real_,
      .groups = "drop")
}
kt <- call_chrom(segs) %>% dplyr::select(-theta)

both_p <- function(n) 1 - 2 * pbinom(1, n, 0.5)     # P(>= 2 reads of each allele | het, n reads, balanced)
bar <- ac[DP >= 10 & OTH == 0, .(CHROM = as.integer(sub("^chr", "", CHROM)),
                                 obs = REF_N >= 2 & ALT_N >= 2, exp = 2 * AF * (1 - AF) * both_p(DP))]
bar <- bar[!is.na(CHROM), .(sites_dp10 = .N, both_allele_ratio = round(sum(obs) / sum(exp), 3)), by = CHROM]
if (opt$allele_source == "vcf") bar[, both_allele_ratio := NA_real_]   # needs every panel site, not just called variants
kt <- kt %>% left_join(bar, by = "CHROM") %>%
  mutate(both_allele_norm = round(both_allele_ratio / median(both_allele_ratio[call == "neu"], na.rm = TRUE), 3),
         copy_change = unname(delta[call])) %>%
  arrange(CHROM)
fwrite(kt, file.path(sd, "karyotype.tsv"), sep = "\t")

n_chrom <- 46 + sum(kt$copy_change, na.rm = TRUE)
group <- cut(n_chrom, c(-Inf, 23, 30, 39, 43, 45, 46, 50, 67, 79, Inf),
             c("<24", "near-haploid (24-30)", "low hypodiploid (31-39)", "high hypodiploid (40-43)",
               "near-diploid hypodiploid (44-45)", "diploid (46)", "low hyperdiploid (47-50)",
               "high hyperdiploid (51-67)", "near-triploid (68-79)", ">79"))
gains <- with(kt, c(rep(CHROM[call == "amp"], 1), rep(CHROM[call == "bamp"], 2)))
fmt <- function(x, sign) if (length(x)) paste0(sign, sort(x), collapse = ",") else ""
n_cnloh <- sum(kt$call == "loh")
flags <- c(
  if (n_cnloh >= opt$cnloh_flag) glue("possible_doubled_hypodiploid({n_cnloh} whole-chromosome CN-LOH)"),
  if (nrow(df_allele) < 5000) glue("low_allele_data({nrow(df_allele)} het sites)"),
  if (any(kt$call_frac < 0.8)) "partial_chromosome_events"
)
summary <- data.table(
  sample = sample,
  est_chromosomes = n_chrom,
  ploidy_group = as.character(group),
  gains = fmt(gains, "+"),
  losses = fmt(kt$CHROM[kt$call %in% c("del", "bdel")], "-"),
  cn_loh = fmt(kt$CHROM[kt$call == "loh"], "UPD"),
  het_sites = nrow(df_allele),
  reference = ref$source,
  mode = "numbat expression only",
  flags = paste(flags, collapse = ";"),
  note = "autosomes only; sex chromosomes assumed 2 (check X with RNAseqCNV)"
)
fwrite(summary, file.path(sd, "summary.tsv"), sep = "\t")
print(summary[, .(sample, est_chromosomes, ploidy_group, gains, losses, cn_loh, het_sites, flags)])
