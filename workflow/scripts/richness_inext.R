#!/usr/bin/env Rscript
# =============================================================================
# richness_inext.R: NEW rarefaction / extrapolation of MOTU richness (iNEXT,
# Hsieh, Ma & Chao 2016).
#
# For every site (and "All") with at least 2 MOTUs:
#   - abundance-based curves (specimens = individuals), q = 0 (richness),
#     extrapolated to endpoint_factor x the observed number of specimens
#   - sample coverage of the observed data
#   - asymptotic richness (Chao1) with 95% confidence interval
#   - incidence-based curves and Chao2 (sampling events = sampling units),
#     for sites with at least 3 sampling events
# Writes a table of estimates, a PDF of curves and metrics.
# =============================================================================

suppressPackageStartupMessages({
  library(iNEXT)
  library(ggplot2)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--abundance"), make_option("--incidence"),
  make_option("--nboot", type = "integer", default = 999),
  make_option("--endpoint_factor", type = "double", default = 2),
  make_option("--out_estimates"), make_option("--out_pdf"), make_option("--out_metrics")
)))

# ---- Abundance data --------------------------------------------------------------
ab <- read.delim(opt$abundance, row.names = 1, check.names = FALSE)
# Keep sites with at least 2 MOTUs (iNEXT needs some data to extrapolate).
ab <- ab[, colSums(ab > 0) >= 2, drop = FALSE]
ab_list <- lapply(ab, function(x) x[x > 0])

out_rows <- list()
metrics <- c()
pdf(opt$out_pdf, width = 7, height = 5)

if (length(ab_list) > 0) {
  endpoint <- ceiling(opt$endpoint_factor * max(sapply(ab_list, sum)))
  res <- iNEXT(ab_list, q = 0, datatype = "abundance", endpoint = endpoint, nboot = opt$nboot)
  print(ggiNEXT(res, type = 1) + ggtitle("MOTU richness vs. number of specimens"))
  print(ggiNEXT(res, type = 2) + ggtitle("Sample coverage vs. number of specimens"))
  print(ggiNEXT(res, type = 3) + ggtitle("MOTU richness vs. sample coverage"))

  # Asymptotic richness (Chao1) and coverage of the observed sample.
  for (site in names(ab_list)) {
    x <- ab_list[[site]]
    est <- ChaoRichness(x, datatype = "abundance")
    cov <- DataInfo(x, datatype = "abundance")
    sc <- if ("SC" %in% names(cov)) cov$SC else NA
    out_rows[[length(out_rows) + 1]] <- data.frame(
      site = site, data = "abundance", n = sum(x), observed = est$Observed,
      estimator = round(est$Estimator, 1), ci_low = round(est$`95% Lower`, 1),
      ci_high = round(est$`95% Upper`, 1), coverage = round(sc, 3))
    s <- gsub("[^A-Za-z0-9]", "_", site)
    metrics[paste0("richness_", s, "_specimens")] <- sum(x)
    metrics[paste0("richness_", s, "_observed")] <- est$Observed
    metrics[paste0("richness_", s, "_chao1")] <- round(est$Estimator, 1)
    metrics[paste0("richness_", s, "_chao1_ci")] <- paste0(round(est$`95% Lower`, 1), "-", round(est$`95% Upper`, 1))
    metrics[paste0("richness_", s, "_coverage")] <- round(sc, 3)
  }
}

# ---- Incidence data -----------------------------------------------------------------
inc <- read.delim(opt$incidence, row.names = 1, check.names = FALSE)
# First row = number of sampling events; need >= 3 events and >= 2 MOTUs.
inc <- inc[, inc[1, ] >= 3 & colSums(inc[-1, , drop = FALSE] > 0) >= 2, drop = FALSE]
if (ncol(inc) > 0) {
  inc_list <- lapply(inc, function(x) c(x[1], x[-1][x[-1] > 0]))
  res_i <- iNEXT(inc_list, q = 0, datatype = "incidence_freq", nboot = opt$nboot)
  print(ggiNEXT(res_i, type = 1) + ggtitle("MOTU richness vs. number of sampling events"))
  for (site in names(inc_list)) {
    est <- ChaoRichness(inc_list[[site]], datatype = "incidence_freq")
    out_rows[[length(out_rows) + 1]] <- data.frame(
      site = site, data = "incidence", n = inc_list[[site]][1], observed = est$Observed,
      estimator = round(est$Estimator, 1), ci_low = round(est$`95% Lower`, 1),
      ci_high = round(est$`95% Upper`, 1), coverage = NA)
    s <- gsub("[^A-Za-z0-9]", "_", site)
    metrics[paste0("richness_", s, "_chao2")] <- round(est$Estimator, 1)
    metrics[paste0("richness_", s, "_sampling_events")] <- inc_list[[site]][1]
  }
}
dev.off()

est_table <- if (length(out_rows)) do.call(rbind, out_rows) else data.frame()
write.table(est_table, opt$out_estimates, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(paste0(names(metrics), "\t", metrics), opt$out_metrics)
