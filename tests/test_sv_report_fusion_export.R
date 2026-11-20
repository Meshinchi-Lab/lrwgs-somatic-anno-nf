#!/usr/bin/env Rscript
# tests/test_sv_report_fusion_export.R
#
# Runs the `sv-proteinpaint-build` chunk of bin/sv_report_index.qmd on a pass-2 table
# and checks the svfusion table and the Gene Fusion app lines. Run from the repository
# root in the report image:
#   docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w \
#     ghcr.io/meshinchi-lab/sv_report_index:2026.11 \
#     Rscript tests/test_sv_report_fusion_export.R [cohort_sv_pass2.tsv]
#
# tests/fixtures/sv_report/fusion_events.tsv (two samples):
#   S1 sv_del_mm      DEL, both genes on "-": the END gene is 5' and comes first
#   S1 sv_inv_2/_1    BND pair, both ALTs t]p] (same orientation), listed _2 first:
#                     the "+" gene kept LEFT is 5' (GENEP), whichever mate is read
#   S2 sv_one_2       only mate present, BCL11B in a gene, TLX3 side intergenic:
#                     the in-gene breakend first
#   S2 sv_anti_2/_1   BND pair joining two genes antisense, listed _2 first: the _1
#                     record's own breakend first; no Gene Fusion line
#   S2 sv_none_2      both breakends intergenic: the record's own breakend first
# With a cohort_sv_pass2.tsv argument (the report writes pass2_clean to it) the
# shuffle and agreement checks also run on that table.

qmd   <- readLines("bin/sv_report_index.qmd")
start <- grep("^```\\{r sv-proteinpaint-build\\}", qmd)
stopifnot(length(start) == 1)
end   <- start + which(qmd[(start + 1):length(qmd)] == "```")[1]
chunk <- parse(text = qmd[(start + 1):(end - 1)])

read_pass2 <- function(path) readr::read_tsv(path, col_types = readr::cols(.default = "c"),
                                             progress = FALSE)

run_chunk <- function(tbl) {
  dir <- tempfile("svfusion_"); dir.create(dir)
  old <- setwd(dir); on.exit(setwd(old))
  env <- new.env(parent = globalenv())
  env$pass2_clean <- tbl
  eval(chunk, envir = env)
  rd  <- function(f) { p <- file.path(dir, "protein_paint_out", f); if (file.exists(p)) readLines(p) else character(0) }
  list(svfusion = rd("sv_tier1_2.svfusion.txt"), genefusion = rd("sv_tier1_2.genefusion.txt"),
       table = env$pp_sv_fusion)
}

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop("FAIL: ", msg, call. = FALSE)
  cat("ok -", msg, "\n")
}

# Each Gene Fusion line names the same breakpoint first as its svfusion row, except
# the second junction of an INV (the svfusion row follows the first junction).
agree <- function(out) {
  t  <- out$table
  sa <- paste(t$chr_a, t$position_a); sb <- paste(t$chr_b, t$position_b)
  all(vapply(strsplit(out$genefusion, "::", fixed = TRUE), function(p) {
    a  <- strsplit(p[1], ",")[[1]]; b <- strsplit(p[2], ",")[[1]]
    ka <- paste(a[2], as.numeric(a[3]) - 1); kb <- paste(b[2], as.numeric(b[3]) - 1)
    hit <- which((sa == ka & sb == kb) | (sa == kb & sb == ka))
    length(hit) >= 1 && all(t$SVTYPE[hit] == "INV" | sa[hit] == ka)
  }, logical(1)))
}

shuffled_identical <- function(tbl, out, label) {
  set.seed(20261120)
  for (k in 1:5) {
    s <- run_chunk(tbl[sample(nrow(tbl)), ])
    check(identical(s$svfusion, out$svfusion) && identical(s$genefusion, out$genefusion),
          sprintf("%s: shuffled rows (%d) give identical svfusion and Gene Fusion files", label, k))
  }
}

# ── fixture ────────────────────────────────────────────────────────────────
fx  <- read_pass2("tests/fixtures/sv_report/fusion_events.tsv")
out <- run_chunk(fx)
t   <- out$table
row <- function(s, chr) t[t$Sample == s & t$chr_a == chr, ]

check(nrow(t) == 5, "5 svfusion events: each BND pair collapsed to one row")
r <- row("S1", "chr1")
check(nrow(r) == 1 && r$gene_a == "GENER" && r$position_a == 4999 && r$gene_b == "GENEL" && r$position_b == 999,
      "DEL of two '-' strand genes: the END gene (5') first")
r <- row("S1", "chr5")
check(nrow(r) == 1 && r$gene_a == "GENEP" && r$strand_a == "+" && r$gene_b == "GENEM" && r$chr_b == "chr14",
      "same-orientation BND pair: the '+' gene kept LEFT (5') first, not the first row's")
r <- row("S2", "chr14")
check(nrow(r) == 1 && r$gene_a == "BCL11B" && r$gene_b == "" && r$chr_b == "chr5",
      "one intergenic breakend: the in-gene breakend first")
r <- row("S2", "chr3")
check(nrow(r) == 1 && r$gene_a == "GENEA" && r$chr_b == "chr4",
      "antisense two-gene BND pair: the _1 record's own breakend first")
r <- row("S2", "chr7")
check(nrow(r) == 1 && r$gene_a == "" && r$gene_b == "" && r$chr_b == "chr8",
      "both breakends intergenic: the record's own breakend first")
check(identical(out$genefusion, c(
        "GENER,chr1,5000,-::GENEL,chr1,1000,-",
        "GENEP,chr5,200,+::GENEM,chr14,100,-",
        "BCL11B,chr14,99232731,-::TLX3,chr5,171319244,+")),
      "Gene Fusion lines: 5' gene first; antisense and two-intergenic events left out")
check(agree(out), "fixture: svfusion and Gene Fusion name the same breakpoint first")
shuffled_identical(fx, out, "fixture")

# ── optional: a real cohort_sv_pass2.tsv ─────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  real <- read_pass2(args[1])
  rout <- run_chunk(real)
  check(nrow(rout$table) > 0, sprintf("%s: %d svfusion events", args[1], nrow(rout$table)))
  check(agree(rout), sprintf("%s: svfusion and Gene Fusion name the same breakpoint first", args[1]))
  shuffled_identical(real, rout, args[1])
}
cat("all checks passed\n")
