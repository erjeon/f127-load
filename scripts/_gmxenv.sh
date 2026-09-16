# Sourced by 3_equilibrate.sh, 4_production.sh and 5_analyze.sh.
# Works out which gmx to call, how many threads it may use, and whether that
# gmx can talk to a GPU at all.
#
# The gmx already on PATH is used, and a GMXRC is sourced only when there is
# none. Which one is in use gets printed.
set +u
command -v gmx >/dev/null 2>&1 || source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true
set -uo pipefail
GMX=${GMX:-gmx}
# PYTHON= picks the interpreter, as in 0_check.sh.
PY=${PYTHON:-python3}
GMX_VERSION_TEXT="$($GMX --version 2>/dev/null || true)"

# A GROMACS built without GPU support stops on -nb gpu, so the flags follow
# what gmx --version reports.
if [[ -z ${GMX_GPU:-} ]]; then
  if grep -qiE '^GPU support: *(CUDA|OpenCL|SYCL|HIP)' <<<"$GMX_VERSION_TEXT"
  then GMX_GPU=yes; else GMX_GPU=no; fi
fi

# Thread count from the machine, capped at 16.
if [[ -z ${NTOMP:-} ]]; then
  NTOMP=$( (command -v nproc >/dev/null 2>&1 && nproc) \
           || sysctl -n hw.ncpu 2>/dev/null || echo 8 )
  (( NTOMP > 16 )) && NTOMP=16
fi

# C-rescale needs GROMACS 2021. On an older gmx the sed that sets nsteps also
# swaps the barostat to Parrinello-Rahman.
GMX_MAJOR=$(grep -m1 -oE '[0-9]{4}' <<<"$(grep -m1 -i 'GROMACS version' <<<"$GMX_VERSION_TEXT")")
if [[ -n ${GMX_MAJOR:-} && $GMX_MAJOR -lt 2021 ]]; then
  BAROSTAT_SED='s/^pcoupl  *=.*C-rescale/pcoupl                  = Parrinello-Rahman/'
else
  BAROSTAT_SED=''
fi

# $1 is em or md. Minimisation gains nothing from PME on the GPU.
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
