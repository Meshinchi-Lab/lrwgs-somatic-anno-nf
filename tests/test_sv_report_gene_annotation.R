#!/usr/bin/env Rscript
# tests/test_sv_report_gene_annotation.R
#
# Runs the EnsDb breakend-gene block of bin/sv_report_index.qmd (chunk `pass2-clean`,
# from "# Optional EnsDb-based BND breakend gene annotation" up to "# has_mate") on a
# two-sample table and checks that each record gets the genes at ITS OWN breakpoints.
# SV caller IDs are per-sample counters, so the same ID names different events in
# different samples (seen in a 52-sample cohort: 44 IDs shared at different loci).
# Run from the repository root in the report image, with the EnsDb mounted:
#   docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w \
#     -v "$(cd data/reference && pwd -P)":/ref \
#     ghcr.io/meshinchi-lab/sv_report_index:2026.11 \
#     Rscript tests/test_sv_report_gene_annotation.R /ref/EnsDb.Hsapiens.v110.sqlite

suppressPackageStartupMessages({ library(dplyr); library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1, file.exists(args[1]))

qmd   <- readLines("bin/sv_report_index.qmd")
start <- grep("^# Optional EnsDb-based BND breakend gene annotation", qmd)
end   <- grep("^# has_mate", qmd)
stopifnot(length(start) == 1, length(end) == 1, end > start)
block <- parse(text = qmd[start:(end - 1)])

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop("FAIL: ", msg, call. = FALSE)
  cat("ok -", msg, "\n")
}

# The same ID in two samples at different loci (genes from EnsDb v110):
#   S1  chr21:34816131 (RUNX1) | chr2:99727265 (AFF3)
#   S2  chr14:99232731 (BCL11B) | chr5:171319244 (intergenic, near TLX3)
pass2_clean <- tibble::tibble(
  SAMPLE     = c("S1", "S2"),
  ID         = c("severus_BND1_1", "severus_BND1_1"),
  SVTYPE     = c("BND", "BND"),
  breakpoint = c("chr21:34816131|chr2:99727265", "chr14:99232731|chr5:171319244")
)
env <- new.env(parent = globalenv())
env$pass2_clean <- pass2_clean
env$params <- list(ensembl_db = args[1])
env$install_bioc_if_missing <- function(pkg) invisible(NULL)
eval(block, envir = env)
out <- env$pass2_clean

check(nrow(out) == 2, "one row per input record")
s1 <- out[out$SAMPLE == "S1", ]; s2 <- out[out$SAMPLE == "S2", ]
check(s1$gene_name_left == "RUNX1" && s1$gene_name_right == "AFF3",
      sprintf("S1 genes from its own breakpoints (got %s / %s)", s1$gene_name_left, s1$gene_name_right))
check(s2$gene_name_left == "BCL11B" && s2$gene_source_right == "nearest",
      sprintf("S2 genes from its own breakpoints (got %s / %s)", s2$gene_name_left, s2$gene_source_right))
cat("all checks passed\n")
