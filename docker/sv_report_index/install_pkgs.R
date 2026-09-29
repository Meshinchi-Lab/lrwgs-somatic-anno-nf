#!/usr/bin/env Rscript
# Resilient package install for the SV_REPORT_INDEX image.
#
# WHY NOT `install2.r` / a bare `BiocManager::install()`:
# neither retries, and neither falls back to a second repository. The build
# failure on 2026-09-17 was NOT an outage —
#
#   status was 'Couldn't resolve host name'   (rspm-sync.rstudio.com)
#
# Measured the same day: `packagemanager.posit.co` resolved on 4 of 5 attempts
# from the default resolver, and resolved reliably via 1.1.1.1 / 8.8.8.8 with
# the P3M PACKAGES index returning HTTP 200. So the service was up and the DNS
# was flaky. A single-shot install turns a transient lookup failure into a
# failed build.
#
# Strategy: try each repo in order, retrying with backoff, then verify. P3M
# first because it serves prebuilt Linux binaries (minutes, not a ~30-minute
# source compile of the Bioconductor annotation tree); CRAN/Bioconductor
# mirrors second because they stay reachable when P3M DNS wobbles.
#
# Usage:  Rscript install_pkgs.R cran DT plotly
#         Rscript install_pkgs.R bioc ensembldb plyranges

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: install_pkgs.R <cran|bioc> <pkg> [pkg...]")
mode <- args[1]
pkgs <- args[-1]

ATTEMPTS <- 4L

# The image's own repos first (P3M in rocker/verse), then public mirrors.
image_repos <- local({
  r <- getOption("repos")
  unname(r[!is.na(r) & nzchar(r) & r != "@CRAN@"])
})
cran_chain <- unique(c(image_repos,
                       "https://cloud.r-project.org",
                       "https://cran.rstudio.com"))

with_retry <- function(label, fn) {
  for (i in seq_len(ATTEMPTS)) {
    ok <- tryCatch({ fn(); TRUE },
                   error   = function(e) { message("  [", label, "] attempt ", i,
                                                   " error: ", conditionMessage(e)); FALSE },
                   warning = function(w) { message("  [", label, "] attempt ", i,
                                                   " warning: ", conditionMessage(w)); FALSE })
    if (ok) return(TRUE)
    if (i < ATTEMPTS) Sys.sleep(5 * i)   # linear backoff: 5s, 10s, 15s
  }
  FALSE
}

missing_of <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

if (mode == "cran") {
  for (repo in cran_chain) {
    todo <- missing_of(pkgs)
    if (!length(todo)) break
    message("== CRAN-type repo: ", repo, "  (need: ", paste(todo, collapse = ", "), ")")
    with_retry(repo, function() install.packages(todo, repos = repo, Ncpus = max(1L, parallel::detectCores())))
  }
} else if (mode == "bioc") {
  # BiocManager itself may be absent; fetch it through the same chain.
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    for (repo in cran_chain) {
      if (requireNamespace("BiocManager", quietly = TRUE)) break
      message("== BiocManager from: ", repo)
      with_retry(repo, function() install.packages("BiocManager", repos = repo))
    }
  }
  if (!requireNamespace("BiocManager", quietly = TRUE))
    stop("BiocManager could not be installed from any of: ",
         paste(cran_chain, collapse = ", "))
  # BiocManager resolves the Bioconductor release matching this R version.
  todo <- missing_of(pkgs)
  if (length(todo)) {
    message("== Bioconductor (release ", as.character(BiocManager::version()),
            "): ", paste(todo, collapse = ", "))
    with_retry("bioc", function()
      BiocManager::install(todo, ask = FALSE, update = FALSE,
                           Ncpus = max(1L, parallel::detectCores())))
  }
} else stop("unknown mode: ", mode)

still <- missing_of(pkgs)
if (length(still)) {
  stop("FAILED to install: ", paste(still, collapse = ", "), "\n",
       "Repos tried: ", paste(cran_chain, collapse = ", "), "\n",
       "If this is DNS flakiness rather than a real outage, give the Docker\n",
       "daemon explicit resolvers and rebuild:\n",
       "  Docker Desktop > Settings > Docker Engine > add\n",
       '      "dns": ["1.1.1.1", "8.8.8.8"]\n',
       "  then Apply & Restart.")
}
for (p in pkgs) message("  ok  ", p, "  ", as.character(packageVersion(p)))
