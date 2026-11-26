#!/usr/bin/env Rscript
# Gene counts from a BAM, htseq-count style, so the pipeline needs only a BAM: union mode (a read counts for a gene
# when every overlapping exon belongs to that one gene), unique alignments only (NH = 1), MAPQ >= 10, no secondary /
# supplementary / QC-fail; duplicates kept (as htseq-count). Paired-end: one count per fragment, from read 1.
# Annotation: Ensembl GRCh38 release 102, the gene set behind RNAseqCNV's diploid reference (all 60,616 IDs of the
# clinical HTSeq files); contig names are matched to the BAM (1 / chr1, MT / chrM).
# Strandedness: --strand auto counts chr22 both ways and picks forward / reverse when one orientation has > 4x the
# other, else unstranded.
# Output: HTSeq format (ENSG <tab> count, no header), every annotated gene.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicAlignments)
  library(Rsamtools)
  library(rtracklayer)
  library(parallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bam", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--exons_rds", type = "character", help = "cache of exons by gene (built from --gtf when missing)"),
  make_option("--out", type = "character"),
  make_option("--threads", type = "integer", default = 4),
  make_option("--strand", type = "character", default = "auto", help = "auto, no, forward or reverse")
)))

if (!is.null(opt$exons_rds) && file.exists(opt$exons_rds)) {
  ex <- readRDS(opt$exons_rds)
} else {
  message("reading exons from ", opt$gtf)
  g <- import(opt$gtf, feature.type = "exon", colnames = "gene_id")
  ex <- reduce(split(g, g$gene_id))
  if (!is.null(opt$exons_rds)) saveRDS(ex, opt$exons_rds)
}

bf <- BamFile(opt$bam)
len <- seqlengths(bf)
# match annotation contig names to the BAM
chr_style <- any(startsWith(names(len), "chr"))
sl <- seqlevels(ex)
new <- if (chr_style) ifelse(startsWith(sl, "chr"), sl, paste0("chr", sub("^MT$", "M", sl))) else sub("^chr", "", sub("^chrM$", "MT", sl))
seqlevels(ex) <- new
gene_chrom <- vapply(seqnames(ex), function(s) as.character(s[1]), "")
chroms <- intersect(names(len), unique(gene_chrom))

# paired-end?
open(bf); flags <- scanBam(BamFile(opt$bam, yieldSize = 2000), param = ScanBamParam(what = "flag"))[[1]]$flag; close(bf)
paired <- any(bitwAnd(flags, 1L) > 0)

read_chrom <- function(ch) {
  p <- ScanBamParam(which = GRanges(ch, IRanges(1, len[[ch]])), mapqFilter = 10, tag = "NH",
                    flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                                       isSupplementaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                                       isFirstMateRead = if (paired) TRUE else NA))
  r <- readGAlignments(opt$bam, param = p)
  r[is.na(mcols(r)$NH) | mcols(r)$NH == 1]
}
count_chrom <- function(ch, strand) {
  r <- read_chrom(ch)
  if (strand == "reverse") r <- invertStrand(r)
  e <- ex[gene_chrom == ch]
  assay(summarizeOverlaps(e, r, mode = "Union", ignore.strand = strand == "no", inter.feature = TRUE))[, 1]
}

strand <- opt$strand
if (strand == "auto") {
  test <- if ("chr22" %in% chroms) "chr22" else if ("22" %in% chroms) "22" else chroms[1]
  fw <- sum(count_chrom(test, "forward")); rv <- sum(count_chrom(test, "reverse"))
  strand <- if (fw > 4 * rv) "forward" else if (rv > 4 * fw) "reverse" else "no"
  message(sprintf("strandedness on %s: forward %d, reverse %d -> %s", test, fw, rv, strand))
}
message(sprintf("counting %d contigs, %s, %s, %d threads", length(chroms), if (paired) "paired (read 1)" else "single-end",
                strand, opt$threads))
counts <- unlist(mclapply(chroms, count_chrom, strand = strand, mc.cores = opt$threads, mc.preschedule = FALSE),
                 use.names = TRUE)
out <- data.table(gene = names(ex), count = 0L)
out[match(names(counts), gene), count := as.integer(counts)]
fwrite(out[order(gene)], opt$out, sep = "\t", col.names = FALSE)
message(sprintf("wrote %s: %d genes, %d reads assigned (%s)", opt$out, nrow(out), sum(out$count), strand))
