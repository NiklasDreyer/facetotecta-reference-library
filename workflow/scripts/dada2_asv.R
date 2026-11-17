#!/usr/bin/env Rscript
# =============================================================================
# dada2_asv.R: amplicon sequence variants (ASVs) with DADA2 (Callahan et al. 2016).
#
# Follows the DADA2 tutorial (https://benjjneb.github.io/dada2/):
#   1. filterAndTrim: remove reads with Ns, expected errors > maxEE,
#      truncate at the first base with quality <= truncQ, remove PhiX
#   2. learnErrors: estimate the error model from the data
#   3. dada: denoise each sample
#   4. mergePairs: merge forward and reverse reads (or forward reads only)
#   5. makeSequenceTable + removeBimeraDenovo: remove chimeras
# Outputs: ASV FASTA (IDs <prefix><n>, ordered by total abundance), an ASV x
# sample count table, a table tracking reads through each step, and metrics.
#
# Used both for our Illumina libraries and for public datasets.
# File lists are passed comma-separated. "--r2 NONE" = single-end reads.
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
})

# ---- Command-line options ----------------------------------------------------
opt_list <- list(
  make_option("--r1", type = "character", help = "comma-separated R1 files"),
  make_option("--r2", type = "character", default = "NONE", help = "comma-separated R2 files, or NONE"),
  make_option("--samples", type = "character", help = "comma-separated sample names (same order)"),
  make_option("--filtdir", type = "character", help = "folder for filtered reads"),
  make_option("--merge", type = "character", default = "TRUE", help = "merge read pairs (TRUE/FALSE)"),
  make_option("--trunc_f", type = "integer", default = 0),
  make_option("--trunc_r", type = "integer", default = 0),
  make_option("--maxee_f", type = "double", default = 2),
  make_option("--maxee_r", type = "double", default = 2),
  make_option("--truncq", type = "integer", default = 2),
  make_option("--min_overlap", type = "integer", default = 12),
  make_option("--chimera", type = "character", default = "consensus"),
  make_option("--prefix", type = "character", default = "ASV"),
  make_option("--label", type = "character", default = "dada2"),
  make_option("--threads", type = "integer", default = 1),
  make_option("--out_fasta", type = "character"),
  make_option("--out_table", type = "character"),
  make_option("--out_track", type = "character"),
  make_option("--out_metrics", type = "character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Split the comma-separated lists.
r1 <- strsplit(opt$r1, ",")[[1]]
samples <- strsplit(opt$samples, ",")[[1]]
paired <- opt$r2 != "NONE" && nchar(opt$r2) > 0
r2 <- if (paired) strsplit(opt$r2, ",")[[1]] else NULL
merge <- toupper(opt$merge) == "TRUE" && paired
stopifnot(length(r1) == length(samples))
dir.create(opt$filtdir, recursive = TRUE, showWarnings = FALSE)

# Single-end public runs get an empty placeholder R2 file from run_cutadapt.sh;
# treat empty R2 files as single-end.
if (paired && all(file.size(r2) == 0)) {
  paired <- FALSE
  merge <- FALSE
}

# ---- 1. Filter and trim ---------------------------------------------------------
f1 <- file.path(opt$filtdir, paste0(samples, "_F_filt.fastq.gz"))
names(f1) <- samples
if (paired) {
  f2 <- file.path(opt$filtdir, paste0(samples, "_R_filt.fastq.gz"))
  names(f2) <- samples
  out <- filterAndTrim(r1, f1, r2, f2,
                       truncLen = c(opt$trunc_f, opt$trunc_r),
                       maxN = 0, maxEE = c(opt$maxee_f, opt$maxee_r),
                       truncQ = opt$truncq, rm.phix = TRUE,
                       compress = TRUE, multithread = opt$threads)
} else {
  out <- filterAndTrim(r1, f1, truncLen = opt$trunc_f,
                       maxN = 0, maxEE = opt$maxee_f, truncQ = opt$truncq,
                       rm.phix = TRUE, compress = TRUE, multithread = opt$threads)
}
rownames(out) <- samples

# Samples with no reads left after filtering cannot be denoised: drop them.
keep <- file.exists(f1) & out[, 2] > 0
if (!all(keep)) message("Samples without reads after filtering: ", paste(samples[!keep], collapse = ", "))
f1 <- f1[keep]
if (paired) f2 <- f2[keep]

# ---- 2. Learn error rates -------------------------------------------------------
err1 <- learnErrors(f1, multithread = opt$threads)
if (paired) err2 <- learnErrors(f2, multithread = opt$threads)

# ---- 3. Denoise -------------------------------------------------------------------
dd1 <- dada(f1, err = err1, multithread = opt$threads)
if (paired) dd2 <- dada(f2, err = err2, multithread = opt$threads)

# ---- 4. Merge pairs (or use forward reads) ------------------------------------------
if (merge) {
  merged <- mergePairs(dd1, f1, dd2, f2, minOverlap = opt$min_overlap, verbose = TRUE)
  seqtab <- makeSequenceTable(merged)
} else {
  seqtab <- makeSequenceTable(dd1)
}

# makeSequenceTable needs a list when there is only one sample.
# ---- 5. Remove chimeras --------------------------------------------------------------
seqtab_nochim <- removeBimeraDenovo(seqtab, method = opt$chimera,
                                    multithread = opt$threads, verbose = TRUE)

# ---- Track reads through the steps ----------------------------------------------------
getN <- function(x) sum(getUniques(x))
n_denoised <- if (length(f1) > 1) sapply(dd1, getN) else getN(dd1)
n_merged <- if (merge) { if (length(f1) > 1) sapply(merged, getN) else getN(merged) } else n_denoised
track <- data.frame(sample = samples, input = out[, 1], filtered = out[, 2],
                    denoised = 0, merged = 0, nonchim = 0)
track[keep, "denoised"] <- n_denoised
track[keep, "merged"] <- n_merged
track[keep, "nonchim"] <- rowSums(seqtab_nochim)
write.table(track, opt$out_track, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Write ASVs (ordered by total abundance) and the count table ------------------------
tot <- colSums(seqtab_nochim)
ord <- order(tot, decreasing = TRUE)
seqs <- colnames(seqtab_nochim)[ord]
ids <- paste0(opt$prefix, seq_along(seqs))
writeLines(paste0(">", ids, "\n", seqs), opt$out_fasta)

counts <- t(seqtab_nochim[, ord, drop = FALSE])
rownames(counts) <- ids
counts <- data.frame(asv = ids, counts, check.names = FALSE)
write.table(counts, opt$out_table, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Metrics ------------------------------------------------------------------------------
lab <- opt$label
m <- c(
  setNames(track$input, paste0(lab, "_", track$sample, "_reads_input")),
  setNames(track$filtered, paste0(lab, "_", track$sample, "_reads_filtered")),
  setNames(track$nonchim, paste0(lab, "_", track$sample, "_reads_nonchim")),
  setNames(length(ids), paste0(lab, "_asvs"))
)
writeLines(paste0(names(m), "\t", m), opt$out_metrics)
message("ASVs written: ", length(ids))
