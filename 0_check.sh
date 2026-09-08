#!/usr/bin/env bash
# Detect the local GROMACS and Python environment and write compat.json.
# Every later stage reads that file instead of assuming a version.
#
# Pukyong National University / NCHM Lab
# Eunryul Jeon  <qlsguswjs@pukyong.ac.kr>

set -uo pipefail

OUT="${1:-compat.json}"
fail() { printf '[FAIL] %s\n' "$*" >&2; exit 1; }
warn() { printf '[warn] %s\n'  "$*" >&2; }
ok()   { printf '[ ok ] %s\n'  "$*"; }

# ---------------------------------------------------------------- GROMACS ---
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
ok "GROMACS $GMX_VERSION  ($GMX)"

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
GPU_OK=false
MDRUN_HELP=$("$GMX" mdrun -h 2>&1 || true)
if printf '%s\n' "$MDRUN_HELP" | grep -q -- '-nb'; then
    if command -v nvidia-smi >/dev/null 2>&1; then
        GPU_LIST=$(nvidia-smi -L 2>/dev/null || true)
        if [ -n "$GPU_LIST" ]; then
            GPU_OK=true
            GPU_NAME=$(printf '%s\n' "$GPU_LIST" | head -1 | cut -c1-60)
            ok "GPU detected: $GPU_NAME"
        fi
    fi
fi
[ "$GPU_OK" = true ] || warn "no usable GPU found, runs will use the CPU"

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
if [ -z "$MDA_VERSION" ]; then
    cat >&2 <<'MSG'
[FAIL] MDAnalysis is not installed.  It is required, there is no fallback.

    python3 -m venv .venv
    source .venv/bin/activate
    pip install "MDAnalysis>=2.0"

  Recent Linux and macOS installs refuse pip into the system Python
  (PEP 668).  Use a virtual environment as above, or pip install --user.
MSG
    exit 1
fi
"$PY" - "$MDA_VERSION" <<'PYEOF' || fail "MDAnalysis 2.0 or newer is required, found $MDA_VERSION"
import sys
parts = sys.argv[1].split(".")
sys.exit(0 if int(parts[0]) >= 2 else 1)
PYEOF
ok "MDAnalysis $MDA_VERSION"

NUMPY_VERSION=$("$PY" -c 'import numpy; print(numpy.__version__)' 2>/dev/null) \
    || fail "NumPy is not installed"
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
