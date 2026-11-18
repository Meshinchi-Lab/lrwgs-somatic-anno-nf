#!/usr/bin/env Rscript
#
# st17_svinft.R -- same-type ST17 DEL/DUP matches for pass-2 SV records, in the
# AnnotSV "SVincludedInFt" style but with a coverage floor.
#
# Why: AnnotSV receives the ST17 DEL/DUP BEDs as FtIncludedInSV, so a reference entry
# is attached only when it lies ENTIRELY INSIDE the called SV. A precise long-read
# call that sits slightly inside its published counterpart is therefore missed: for
# example the BCL11B-enhancer tandem duplication called as 3H_PBM68877
# severus_DUP2115 (chr14:98,542,293-98,544,400) against ST17 entries at
# 98,541,878-98,544,474. Annotating the other way round (SV inside the reference)
# alone would over-call, because many ST17 entries span megabases; the coverage floor
# keeps an entry only when the SV covers at least `min_ref_coverage` of its length.
#
# Rule, per DEL or DUP record: attach every ST17 entry of the SAME type whose overlap
# with the SV covers >= min_ref_coverage of the entry's length. This includes entries
# inside the SV (coverage 1, already reported by AnnotSV; the report merge de-duplicates)
# and entries that contain or partly overlap the SV.
#
# Output: the input pass-2 table with every original column written back byte for byte,
# plus two columns:
#   ST17_DEL_SVinFt, ST17_DUP_SVinFt   "name [score]" entries joined by "|", ordered by
#                                      score (desc) then name; empty when nothing matches
#                                      or the record is not of that type.
# Names are encoded as AnnotSV writes them into the pass-2 table (";" -> "|", " " -> "_"),
# so the report can merge them with AnnotSV's name_ST17_DEL / name_ST17_DUP cells.
#
# The EnsDb supplies the chromosome names and lengths (Ensembl style: 1..22, X, Y):
# records or entries on other contigs, or past a chromosome end, are not matched.
#
# Usage: st17_svinft.R <in.tsv> <out.tsv> <st17_del.bed> <st17_dup.bed> <ensdb.sqlite> <min_ref_coverage>

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(plyranges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("usage: st17_svinft.R <in.tsv> <out.tsv> <st17_del.bed> <st17_dup.bed> <ensdb.sqlite> <min_ref_coverage>")
}
in_tsv   <- args[1]
out_tsv  <- args[2]
bed_path <- c(DEL = args[3], DUP = args[4])
ensdb    <- args[5]
min_cov  <- suppressWarnings(as.numeric(args[6]))
if (is.na(min_cov) || min_cov <= 0 || min_cov > 1) {
  stop("min_ref_coverage must be a number in (0, 1], got '", args[6], "'")
}
# Never write over the input: in a Nextflow task the input is a symlink into the upstream
# process's output, so writing to it would overwrite that file.
if (file.exists(out_tsv) &&
    normalizePath(out_tsv, mustWork = FALSE) == normalizePath(in_tsv, mustWork = FALSE)) {
  stop("output path resolves to the input file; refusing to overwrite: ", in_tsv)
}

norm_chr <- function(x) sub("^chr", "", x)
# AnnotSV's encoding of a user-BED name inside the pass-2 table.
encode_name <- function(x) gsub(" ", "_", gsub(";", "|", x, fixed = TRUE), fixed = TRUE)

# ---- chromosomes and lengths from the EnsDb ------------------------------------------
con <- DBI::dbConnect(RSQLite::SQLite(), ensdb, flags = RSQLite::SQLITE_RO)
chrom <- DBI::dbGetQuery(con, "SELECT seq_name, seq_length FROM chromosome")
DBI::dbDisconnect(con)
chrom <- chrom[chrom$seq_name %in% c(as.character(1:22), "X", "Y"), , drop = FALSE]
if (nrow(chrom) == 0) stop("no primary chromosomes (1-22, X, Y) in the EnsDb chromosome table: ", ensdb)
chr_len <- setNames(as.numeric(chrom$seq_length), chrom$seq_name)

# ---- pass-2 table, read as raw lines so original columns are written back unchanged --
lines <- readLines(in_tsv, warn = FALSE)
if (length(lines) == 0) stop("empty input table (no header): ", in_tsv)
header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
need <- c("CHROM", "POS", "END", "SVTYPE")
miss <- setdiff(need, header)
if (length(miss) > 0) stop(in_tsv, " lacks required column(s): ", paste(miss, collapse = ", "))
body   <- lines[-1]
fields <- strsplit(body, "\t", fixed = TRUE)
col_of <- function(col) {
  i <- match(col, header)
  vapply(fields, function(f) if (length(f) >= i) f[i] else NA_character_, character(1))
}
sv <- data.frame(
  row      = seq_along(body),
  seqnames = norm_chr(col_of("CHROM")),
  start    = suppressWarnings(as.numeric(col_of("POS"))),
  end      = suppressWarnings(as.numeric(col_of("END"))),
  svtype   = toupper(col_of("SVTYPE")),
  stringsAsFactors = FALSE
)
sv$len_ok <- sv$seqnames %in% names(chr_len) & !is.na(sv$start) & !is.na(sv$end) &
             sv$start >= 1 & sv$end >= sv$start
sv$len_ok[sv$len_ok] <- sv$end[sv$len_ok] <= chr_len[sv$seqnames[sv$len_ok]]

# ---- ST17 reference entries ----------------------------------------------------------
read_ref <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x) & !startsWith(x, "#")]
  f <- strsplit(x, "\t", fixed = TRUE)
  if (length(f) == 0) return(data.frame())
  if (any(lengths(f) < 5)) stop(path, ": every BED line needs chrom, start, end, name and score")
  d <- data.frame(
    seqnames = norm_chr(vapply(f, `[`, character(1), 1)),
    start    = as.numeric(vapply(f, `[`, character(1), 2)) + 1,   # BED 0-based -> 1-based
    end      = as.numeric(vapply(f, `[`, character(1), 3)),
    name     = encode_name(vapply(f, `[`, character(1), 4)),
    score    = vapply(f, `[`, character(1), 5),
    stringsAsFactors = FALSE
  )
  keep <- d$seqnames %in% names(chr_len) & d$end >= d$start
  d <- d[keep, , drop = FALSE]
  d$ref_width <- d$end - d$start + 1
  d
}

# ---- matching -------------------------------------------------------------------------
out_cols <- list()
for (t in c("DEL", "DUP")) {
  cell <- rep("", length(body))
  ref  <- read_ref(bed_path[[t]])
  q    <- sv[sv$len_ok & sv$svtype == t, c("seqnames", "start", "end", "row"), drop = FALSE]
  n_hit <- 0L
  if (nrow(q) > 0 && nrow(ref) > 0) {
    hits <- suppressWarnings(join_overlap_intersect(as_granges(q), as_granges(ref)))
    if (length(hits) > 0) {
      h <- as.data.frame(hits)
      h <- h[h$width / h$ref_width >= min_cov, , drop = FALSE]
      if (nrow(h) > 0) {
        h$score_num <- suppressWarnings(as.numeric(h$score))
        h <- h[order(h$row, -h$score_num, h$name), , drop = FALSE]
        h <- h[!duplicated(h[, c("row", "name")]), , drop = FALSE]
        agg <- tapply(paste0(h$name, " [", h$score, "]"), h$row, paste, collapse = "|")
        cell[as.integer(names(agg))] <- unname(agg)
        n_hit <- length(agg)
      }
    }
  }
  message(sprintf("st17_svinft: %s: %d records, %d ST17 entries, %d records matched (coverage >= %g)",
                  t, nrow(q), nrow(ref), n_hit, min_cov))
  out_cols[[t]] <- cell
}
n_skip <- sum(sv$svtype %in% c("DEL", "DUP") & !sv$len_ok)
if (n_skip > 0) {
  message(sprintf("st17_svinft: %d DEL/DUP records not matched (contig outside 1-22/X/Y, unparseable or out-of-range coordinates)", n_skip))
}

# ---- write: original lines unchanged, two columns appended ----------------------------
out <- c(paste(lines[1], "ST17_DEL_SVinFt", "ST17_DUP_SVinFt", sep = "\t"),
         if (length(body) > 0) paste(body, out_cols$DEL, out_cols$DUP, sep = "\t") else character(0))
writeLines(out, out_tsv)
