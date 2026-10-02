# shellcheck shell=bash
#
# Sourced by bin/test_run.sh and bin/main_run.sh. Resolves a Nextflow launcher
# pinned to $NXF_PINNED_VER and exports $NXF_CMD for the caller to invoke.
# manifest.nextflowVersion = '!25.10.0' in nextflow.config is the authoritative
# guard and will hard-fail in that situation, but only after JVM startup and
# with a message that does not explain why the export had no effect. Checking
# here fails faster and says what to do about it.

NXF_PINNED_VER=25.10.0
export NXF_VER="$NXF_PINNED_VER"

# NXF_LAUNCHER points at an NXF_VER-aware launcher without changing the
# `nextflow` on PATH, so a Homebrew install can stay in place for other projects.
NXF_CMD="${NXF_LAUNCHER:-nextflow}"

if ! command -v "$NXF_CMD" >/dev/null 2>&1 && [ ! -x "$NXF_CMD" ]; then
    echo "ERROR: Nextflow launcher '$NXF_CMD' not found on PATH and not executable." >&2
    exit 1
fi

_nxf_reported="$("$NXF_CMD" -v 2>/dev/null || true)"

case "$_nxf_reported" in
*"$NXF_PINNED_VER"*) ;;
*)
    cat >&2 <<EOF
ERROR: this pipeline is pinned to Nextflow $NXF_PINNED_VER, for compatibility with the
       upstream Lumos pipeline (enforced by manifest.nextflowVersion = '!$NXF_PINNED_VER').

         launcher : $(command -v "$NXF_CMD" 2>/dev/null || echo "$NXF_CMD")
         reports  : ${_nxf_reported:-<no version reported>}

       This launcher ignored NXF_VER=$NXF_PINNED_VER, so the pinned version would not be
       the one that runs. Install the official launcher, which honours NXF_VER:

         curl -s https://get.nextflow.io | bash
         mkdir -p "\$HOME/.nextflow/launcher" && mv nextflow "\$HOME/.nextflow/launcher/"

       then re-run, either per-invocation:

         NXF_LAUNCHER="\$HOME/.nextflow/launcher/nextflow" bash bin/$(basename "$0")

       or by exporting NXF_LAUNCHER from your shell profile.
EOF
    exit 1
    ;;
esac

unset _nxf_reported
export NXF_CMD
