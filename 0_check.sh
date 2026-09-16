#!/usr/bin/env bash
# Detect the local GROMACS and Python environment and write compat.json.
# Every later stage reads that file instead of assuming a version.
#
# Pukyong National University / NCHM Lab
# Eunryul Jeon  <qlsguswjs@pukyong.ac.kr>

set -uo pipefail

OUT="${1:-compat.json}"
fail() { printf '[failed] %s\n' "$*" >&2; exit 1; }
warn() { printf '[warn] %s\n'  "$*" >&2; }
ok()   { printf '[ ok ] %s\n'  "$*"; }

# ---------------------------------------------------------------- GROMACS ---
# The run scripts source GMXRC before doing anything, so this check has to do the
# same or it reports no GROMACS on a machine where the pipeline would run fine.
# A login shell on a cluster node often has nothing on PATH.
if ! command -v gmx >/dev/null 2>&1 && [ -z "${GMX:-}" ]; then
    for rc in "${GMXRC:-}" /usr/local/gromacs/bin/GMXRC /usr/local/gromacs-*/bin/GMXRC \
              /opt/gromacs/bin/GMXRC "$HOME/gromacs/bin/GMXRC"; do
        [ -n "$rc" ] && [ -f "$rc" ] && { set +u; . "$rc"; set -u; break; }
    done
fi

GMX="${GMX:-}"
if [ -z "$GMX" ]; then
    for candidate in gmx gmx_mpi gmx_d; do
        command -v "$candidate" >/dev/null 2>&1 && { GMX="$candidate"; break; }
    done
fi
[ -n "$GMX" ] || fail "no GROMACS binary on PATH. Source your GMXRC, or set GMX=/path/to/gmx"

# Note: never pipe into an early-exiting filter (grep -q, grep -m1, head) while
# pipefail is set.  The filter closes the pipe, the producer dies on SIGPIPE and
# pipefail reports the whole pipeline as failed.  Capture first, then filter.
GMX_HELP=$("$GMX" --version 2>/dev/null || true)
VERSION_LINE=$(printf '%s\n' "$GMX_HELP" | grep -i 'GROMACS version' | head -1)
GMX_VERSION=$(printf '%s' "$VERSION_LINE" | grep -oE '[0-9]{4}(\.[0-9]+)*' | head -1)
[ -n "$GMX_VERSION" ] || fail "could not parse a version from: $VERSION_LINE"
GMX_MAJOR=${GMX_VERSION%%.*}
ok "GROMACS $GMX_VERSION  ($(command -v "$GMX" 2>/dev/null || printf '%s' "$GMX"))"

# More than one GROMACS on a machine is normal, and the pipeline uses whichever
# is on PATH. Half a day went into a run where 2022.3 was picked up and refused
# the tpr that 2025.4 had written, so the others are listed here rather than
# left to be discovered from a trjconv error.
OTHERS=""
for rc in /usr/local/gromacs*/bin/GMXRC /opt/gromacs*/bin/GMXRC "$HOME"/gromacs*/bin/GMXRC; do
    [ -f "$rc" ] || continue
    other="${rc%/bin/GMXRC}/bin/gmx"
    [ -x "$other" ] || continue
    [ "$(command -v "$GMX" 2>/dev/null)" = "$other" ] && continue
    v=$("$other" --version 2>/dev/null | grep -i 'GROMACS version' | grep -oE '[0-9]{4}(\.[0-9]+)*' | head -1)
    OTHERS="$OTHERS  $other ($v)"
done
if [ -n "$OTHERS" ]; then
    warn "another GROMACS is installed and will not be used:$OTHERS"
    warn "a tpr written by a newer GROMACS cannot be read by an older one. To pick a different build, source its GMXRC before running, or set GMX=/path/to/gmx"
fi

[ "$GMX_MAJOR" -ge 2019 ] || fail "GROMACS $GMX_VERSION is too old. 2019 or newer is required."
[ "$GMX_MAJOR" -ge 2021 ] || warn "GROMACS $GMX_VERSION has no C-rescale barostat. Falling back to Parrinello-Rahman."

# Feature table.  Each entry is the first major version that provides it.
feature() { [ "$GMX_MAJOR" -ge "$1" ] && echo true || echo false; }
HAS_CRESCALE=$(feature 2021)      # pcoupl = C-rescale
HAS_UPDATE_GPU=$(feature 2020)    # mdrun -update gpu
HAS_BONDED_GPU=$(feature 2019)    # mdrun -bonded gpu
NEEDS_XTCOUT=false                # nstxtcout, only below 5.0, never reached here

# mdp keywords that differ between versions
BAROSTAT="C-rescale"; [ "$HAS_CRESCALE" = true ] || BAROSTAT="Parrinello-Rahman"

# ------------------------------------------------------------------- GPU ----
# Two separate things have to be true. A card has to be present, and this gmx
# has to have been built to use it. The Homebrew build on macOS is a normal
# GROMACS with "GPU support: disabled", and asking it for -nb gpu is a hard
# error rather than a fallback, so the run scripts read the same line.
GPU_OK=false
GMX_GPU_BUILD=$(printf '%s\n' "$GMX_HELP" | grep -iE '^GPU support:' | sed 's/.*: *//')
case "$GMX_GPU_BUILD" in
    CUDA*|OpenCL*|SYCL*|HIP*) GPU_BUILT=true ;;
    *)                        GPU_BUILT=false ;;
esac
if [ "$GPU_BUILT" = true ]; then
    GPU_LIST=$(nvidia-smi -L 2>/dev/null || true)
    if [ -n "$GPU_LIST" ]; then
        GPU_OK=true
        ok "GPU: $(printf '%s\n' "$GPU_LIST" | head -1 | cut -c1-60)  (gmx built with $GMX_GPU_BUILD)"
    else
        warn "this gmx was built with $GMX_GPU_BUILD but no card answered. Runs will use the CPU"
    fi
else
    ok "no GPU, runs will use the CPU (this gmx reports GPU support: ${GMX_GPU_BUILD:-disabled})"
fi

# ---------------------------------------------------------------- Python ----
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || fail "no python3 on PATH"
PY_VERSION=$("$PY" -c 'import sys; print("%d.%d"%sys.version_info[:2])')
"$PY" -c 'import sys; sys.exit(0 if sys.version_info>=(3,8) else 1)' \
    || fail "Python $PY_VERSION is too old, 3.8 or newer is required"
ok "Python $PY_VERSION  ($PY)"

MDA_VERSION=$("$PY" - <<'PYEOF' 2>/dev/null
try:
    import MDAnalysis
    print(MDAnalysis.__version__)
except Exception:
    pass
PYEOF
)
# Only the analysis imports MDAnalysis. Building, minimising, equilibrating and
# running do not, so a missing one is not a reason to report the tool unusable:
# it stops four working stages from being tried. It is still loud, because the
# alternative is finding out after the production run.
if [ -z "$MDA_VERSION" ]; then
    MDA_VERSION=none
    cat >&2 <<MSG
[warn] MDAnalysis is not installed, so \`$PY\` can build and run but not analyse.
       Stages 1 to 4 work. Stage 5, scripts/5_analyze.sh, will stop.

           python3 -m venv venv
           source venv/bin/activate
           pip install "MDAnalysis>=2.8" matplotlib

       Recent Linux and macOS refuse pip into the system Python (PEP 668), so
       the virtual environment is not optional. If one is already set up, point
       at it instead:  PYTHON=/path/to/venv/bin/python ./f127 check
MSG
else
    "$PY" - "$MDA_VERSION" <<'PYEOF' || fail "MDAnalysis 2.0 or newer is required, found $MDA_VERSION"
import sys
parts = sys.argv[1].split(".")
sys.exit(0 if int(parts[0]) >= 2 else 1)
PYEOF
    ok "MDAnalysis $MDA_VERSION"
fi

# GROMACS 2025 writes tpx 137 and MDAnalysis only learned to read it in 2.8.
# The analysis step falls back to the .gro beside the .tpr, which costs nothing
# but exact masses, so this is a note and not a failure.
[ "$MDA_VERSION" = none ] || "$PY" - "$MDA_VERSION" <<'PYEOF' || warn \
  "MDAnalysis $MDA_VERSION cannot read the tpr GROMACS 2025 writes (tpx 137). The analysis will read the .gro instead and guess masses from atom names. Upgrade to 2.8 or newer to use the tpr."
import sys
p = [int(x) for x in sys.argv[1].split(".")[:2] if x.isdigit()]
sys.exit(0 if p >= [2, 8] else 1)
PYEOF

NUMPY_VERSION=$("$PY" -c 'import numpy; print(numpy.__version__)' 2>/dev/null) \
    || fail "NumPy is not installed, and the build needs it. pip install numpy"
ok "NumPy $NUMPY_VERSION"

# ------------------------------------------------------------------ write ---
cat > "$OUT" <<JSON
{
  "gmx": "$GMX",
  "gmx_version": "$GMX_VERSION",
  "gmx_major": $GMX_MAJOR,
  "barostat": "$BAROSTAT",
  "has_crescale": $HAS_CRESCALE,
  "has_update_gpu": $HAS_UPDATE_GPU,
  "has_bonded_gpu": $HAS_BONDED_GPU,
  "gpu": $GPU_OK,
  "python": "$PY",
  "python_version": "$PY_VERSION",
  "mdanalysis_version": "$MDA_VERSION",
  "numpy_version": "$NUMPY_VERSION"
}
JSON
ok "wrote $(cd "$(dirname "$OUT")" && pwd)/$(basename "$OUT")"
