# SV_REPORT_INDEX report-rendering image.
#
# WHY THIS EXISTS
# ---------------
# `bin/sv_report_index.qmd` previously installed DT, plotly, ensembldb and
# plyranges at render time. That puts an outbound network dependency on the
# critical path of every report, on a compute node. 
#
# BUILD
# -----
# From the repository root. The two --build-arg values stamp provenance labels;
# --provenance/--sbom attach BuildKit attestations so the public package carries
# a verifiable record of what was built and from which context.
#
#   docker build --check docker/sv_report_index          # Dockerfile lint, no build
#   docker build \
#     --build-arg IMAGE_REVISION="$(git rev-parse HEAD)" \
#     --build-arg IMAGE_CREATED="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
#     --provenance=true --sbom=true \
#     -t ghcr.io/meshinchi-lab/sv_report_index:2026.10 \
#     docker/sv_report_index
#   docker push ghcr.io/meshinchi-lab/sv_report_index:2026.10
#
# The build context is this directory only. Its .dockerignore admits nothing
# but the Dockerfile and install_pkgs.R (it applies only when this directory
# is the context). A repo-root context (`-f docker/sv_report_index/Dockerfile .`)
# is also safe, for a different reason: `COPY install_pkgs.R` then fails, so the
# build cannot complete with data/, results/ or work/ in reach.
#
# PLATFORM: rocker/verse:4.6.1 publishes linux/amd64 ONLY 
#
# USE
# ---
# Point the module at it (see modules/local/sv_report_index/main.nf), or
# override per-run without editing code:
#   nextflow run main.nf \
#     --process.withName:SV_REPORT_INDEX.container 'ghcr.io/meshinchi-lab/sv_report_index:2026.10'
#
# NOTE ON BUILD TIME
# ------------------
# rocker/verse configures Posit Package Manager (P3M), which serves prebuilt
# Linux binaries — the Bioconductor layer takes minutes rather than the ~30+
# it would take compiling ensembldb's dependency tree from source. If your
# build host cannot reach P3M this will fall back to source builds and be slow,
# but it will still succeed.

# Overridable so this is not a hardcoded constant (buildx lint
# FromPlatformFlagConstDisallowed), but defaulted because the base has no other
# option: `docker build --build-arg BASE_PLATFORM=... docker/sv_report_index`
ARG BASE_PLATFORM=linux/amd64
# Pinned by digest, not tag. The `4.6.1` tag is mutable upstream and was
# rebuilt at least twice between 2026-09-23 and 2026-10-07; a tag-only FROM
# therefore does not reproduce. This digest is the manifest list every render
# ran against during validation (2026-09-30 .. 2026-10-07). Bump deliberately:
#   docker manifest inspect rocker/verse:4.6.1   # new list digest
FROM --platform=${BASE_PLATFORM} rocker/verse:4.6.1@sha256:408c531408bb3db00dacb2e98311d1a4671ba7f128688b3649630a6a480a239e

# Fail the build on any error in a pipeline, not just the last command.
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# System libraries. rocker/verse already carries build-essential and most of
# what the tidyverse needs; these cover the XML/curl/SQLite paths that the
# Bioconductor annotation stack pulls in when a binary is unavailable and it
# has to build from source.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- CRAN layer -----------------------------------------------------------
# Deliberately NOT overriding `repos`: the base image already points at P3M
# binaries. Hardcoding cloud.r-project.org here would throw that away — the
# exact mistake that made the runtime install fail behind an institutional
# proxy.
# install2.r does not retry and does not fall back to a second repo. The
# observed failure was INTERMITTENT DNS, not an outage:
#   status was 'Couldn't resolve host name'  (rspm-sync.rstudio.com)
# and packagemanager.posit.co resolved on only some attempts from the default
# resolver while resolving reliably via 1.1.1.1 / 8.8.8.8. So: retry, and fall
# back to CRAN (verified reachable, HTTP 200) if P3M stays unresolvable.
COPY install_pkgs.R /tmp/install_pkgs.R
RUN Rscript /tmp/install_pkgs.R cran DT plotly \
    && rm -rf /tmp/downloaded_packages

# ---- Bioconductor layer ---------------------------------------------------
# BiocManager resolves the Bioconductor release matching this R version, so no
# Bioc version is pinned here — pinning it would risk a mismatch on a base
# image bump. `ask=FALSE, update=FALSE` keeps the build non-interactive and
# stops it from churning the base image's tidyverse.
RUN Rscript /tmp/install_pkgs.R bioc ensembldb plyranges \
    && rm -rf /tmp/downloaded_packages

# ---- Verification ---------------------------------------------------------
# A silently incomplete image is worse than a failed build: the pipeline would
# get all the way to render before dying. Load every package the qmd needs and
# fail the build here if any is missing.
RUN Rscript -e ' \
      pkgs <- c("tidyverse","DT","plotly","jsonlite","htmltools", \
                "ensembldb","plyranges","S4Vectors","GenomeInfoDb","BiocGenerics"); \
      missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; \
      if (length(missing)) stop("MISSING: ", paste(missing, collapse=", ")); \
      cat("all", length(pkgs), "packages present\n"); \
      for (p in pkgs) cat(sprintf("  %-16s %s\n", p, as.character(packageVersion(p))))'

# Quarto ships with rocker/verse; surface its version at build time so an image
# bump that changes Quarto is visible in the build log rather than at render.
RUN quarto --version


ARG IMAGE_SOURCE=https://github.com/meshinchi-lab/lrwgs-somatic-anno-nf
ARG IMAGE_VERSION=2026.10
# Stamped at build time (see BUILD above); "unknown" marks an unstamped local build.
ARG IMAGE_REVISION=unknown
ARG IMAGE_CREATED=unknown
ARG IMAGE_AUTHORS="Jenny L. Smith, Meshinchi Lab"
LABEL org.opencontainers.image.source="${IMAGE_SOURCE}" \
      org.opencontainers.image.url="${IMAGE_SOURCE}" \
      org.opencontainers.image.documentation="${IMAGE_SOURCE}/blob/main/docker/sv_report_index/Dockerfile" \
      org.opencontainers.image.title="sv_report_index" \
      org.opencontainers.image.description="rocker/verse 4.6.1 + DT, plotly, ensembldb, plyranges for the SV_REPORT_INDEX Quarto report" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.revision="${IMAGE_REVISION}" \
      org.opencontainers.image.created="${IMAGE_CREATED}" \
      org.opencontainers.image.base.name="docker.io/rocker/verse:4.6.1" \
      org.opencontainers.image.base.digest="sha256:408c531408bb3db00dacb2e98311d1a4671ba7f128688b3649630a6a480a239e" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.vendor="Meshinchi Lab" \
      org.opencontainers.image.authors="${IMAGE_AUTHORS}"
