#!/usr/bin/env Rscript
# Diploid B-ALL expression reference for numbat, plus the gene list used to restrict SNP sites.
#
# Default source: RNAseqCNV's diploid_standard_hg38 (44 B-ALL samples without large CNVs, raw HTSeq counts).
# --standard_counts: a text file of in-house HTSeq count files (one path per line) from karyotypically normal samples,
# used instead - same lab pipeline as the test samples, so no cross-lab batch effect.
# Each reference sample is ENSG -> symbol, summed by symbol and normalised to proportions; numbat's fit_ref_sse later
# fits each test sample as a non-negative mixture of these columns.
#
# Output: <out_dir>/ball_reference.rds  (list: lambdas = gene x sample matrix, weight_genes = RNAseqCNV weight_table
#         genes as symbols, source), <out_dir>/genes.tsv (numbat gtf_hg38 gene bodies: CHROM start end gene)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(numbat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rnaseqcnv_dir", type = "character"),
  make_option("--standard_counts", type = "character", default = NULL),
  make_option("--out_dir", type = "character")
)))

ens_to_symbol <- function(ensg) {
  suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*", "", ensg), keytype = "ENSEMBL",
                                         column = "SYMBOL", multiVals = "first"))
}

if (is.null(opt$standard_counts)) {
  load(file.path(opt$rnaseqcnv_dir, "data", "diploid_standard_hg38.rda"))
  ref <- as.data.table(diploid_standard_hg38)
  source <- "RNAseqCNV diploid_standard_hg38"
} else {
  files <- readLines(opt$standard_counts)
  ref <- Reduce(function(a, b) merge(a, b, by = "ENSG"), lapply(files, function(f) {
    d <- fread(f, header = FALSE, col.names = c("ENSG", "count"))[!startsWith(ENSG, "__")]
    setnames(d, "count", sub("\\.HTSeq$", "", basename(f)))
  }))
  source <- paste("in-house:", opt$standard_counts)
}
cols <- setdiff(names(ref), "ENSG")
ref[, symbol := ens_to_symbol(ENSG)]
ref <- ref[!is.na(symbol) & symbol != "", lapply(.SD, sum), by = symbol, .SDcols = cols]
lambdas <- as.matrix(ref[, ..cols]); rownames(lambdas) <- ref$symbol
lambdas <- sweep(lambdas, 2, colSums(lambdas), "/")

sysdata <- new.env(); load(file.path(opt$rnaseqcnv_dir, "R", "sysdata.rda"), envir = sysdata)
weight_genes <- unique(na.omit(ens_to_symbol(sysdata$weight_table$ENSG)))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(lambdas = lambdas, weight_genes = weight_genes, source = source),
        file.path(opt$out_dir, "ball_reference.rds"))
fwrite(gtf_hg38[, .(CHROM = paste0("chr", CHROM), start = gene_start, end = gene_end, gene)],
       file.path(opt$out_dir, "genes.tsv"), sep = "\t")
# the same ENSG -> symbol map for run_numbat.R, so runs need no org.Hs.eg.db (absent from the numbat container)
all_ensg <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
fwrite(data.table(ENSG = all_ensg, symbol = ens_to_symbol(all_ensg))[!is.na(symbol)],
       file.path(opt$out_dir, "ens_to_symbol.tsv"), sep = "\t")
message(sprintf("reference: %d samples x %d genes (%s); %d weight genes", ncol(lambdas), nrow(lambdas), source,
                length(weight_genes)))
