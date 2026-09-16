# Sourced by 3_equilibrate.sh, 4_production.sh and 5_analyze.sh.
# Works out which gmx to call, how many threads it may use, and whether that
# gmx can talk to a GPU at all.
#
# Use the gmx that is already on the PATH. Anyone running MD has set one up, and
# choosing a different one for them is not this script's business: a version is
# sometimes pinned on purpose. A GMXRC is sourced only when there is no gmx at
# all. Which one is in use gets printed, because the failure that started this
# was a silent swap from 2025.4 to 2022.3.
set +u
command -v gmx >/dev/null 2>&1 || source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true
set -uo pipefail
GMX=${GMX:-gmx}
# 0_check.sh honours PYTHON= already. The run scripts called python3 directly,
# so pointing the check at a venv and the run at the system Python was possible.
PY=${PYTHON:-python3}
GMX_VERSION_TEXT="$($GMX --version 2>/dev/null || true)"

# A GROMACS built without CUDA stops at "-nb gpu was requested, but the GROMACS
# binary has been built without GPU support". The Homebrew build on macOS is one
# of those, so the flags cannot be hard-coded. gmx reports what it can do.
if [[ -z ${GMX_GPU:-} ]]; then
  if grep -qiE '^GPU support: *(CUDA|OpenCL|SYCL|HIP)' <<<"$GMX_VERSION_TEXT"
  then GMX_GPU=yes; else GMX_GPU=no; fi
fi

# 16 threads used to be hard-coded, which oversubscribes a laptop and leaves a
# 64-core node idle.
if [[ -z ${NTOMP:-} ]]; then
  NTOMP=$( (command -v nproc >/dev/null 2>&1 && nproc) \
           || sysctl -n hw.ncpu 2>/dev/null || echo 8 )
  (( NTOMP > 16 )) && NTOMP=16
fi

# $1 is em or md. Minimisation gains nothing from PME on the GPU.
# C-rescale arrived in GROMACS 2021 and both mdp files ask for it. On an older
# gmx grompp stops with an invalid enum, so the sed that sets nsteps also swaps
# the barostat. Parrinello-Rahman is the documented substitute and 0_check.sh
# already says which one is in use.
GMX_MAJOR=$(grep -m1 -oE '[0-9]{4}' <<<"$(grep -m1 -i 'GROMACS version' <<<"$GMX_VERSION_TEXT")")
if [[ -n ${GMX_MAJOR:-} && $GMX_MAJOR -lt 2021 ]]; then
  BAROSTAT_SED='s/^pcoupl  *=.*C-rescale/pcoupl                  = Parrinello-Rahman/'
else
  BAROSTAT_SED=''
fi

mdrun_opt() {
  if [[ -n ${MDRUN_OPT:-} ]]; then echo "$MDRUN_OPT"; return; fi
  local o="-ntmpi 1 -ntomp $NTOMP"
  if [[ $GMX_GPU == yes ]]; then
    [[ ${1:-md} == em ]] && o="$o -nb gpu" || o="$o -nb gpu -bonded gpu -pme gpu"
  fi
  echo "$o"
}

gmx_line() {
  echo "  using $(grep -m1 -o 'GROMACS.*' <<<"$GMX_VERSION_TEXT" || echo 'gmx (version unknown)')"
}

gmx_banner() {
  gmx_line
  if [[ $GMX_GPU == yes ]]
  then echo "  $NTOMP threads, work offloaded to the GPU"
  else echo "  $NTOMP threads, on the CPU (this gmx was built without GPU support)"; fi
}
