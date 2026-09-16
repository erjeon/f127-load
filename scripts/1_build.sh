#!/usr/bin/env bash
# Build, and optionally run, the system described by a config from `f127 new`.
#
#   bash scripts/1_build.sh system.cfg            build only, stops at the tpr
#   bash scripts/1_build.sh --run system.cfg 100  build, equilibrate, run 100 ns
#
# route "both" builds two directories from one config so the pair differs only
# in where the solute started, which is the comparison the paper rests on.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=${PYTHON:-python3}
RUN=0; [[ ${1:-} == --run ]] && { RUN=1; shift; }
# The equilibration length was reachable only through NPT_NS, which is not
# somewhere anyone looks. It is an argument now, and the variable still works.
NPT_NS=${NPT_NS:-2}
ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    --npt) NPT_NS=$2; shift 2 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"
NS=${2:-100}
# Naming a file that is not there used to fall through to system.cfg and then to
# system.json in the current directory, so a mistyped path quietly built whatever
# config happened to be lying around. The two defaults apply only when no file
# was named at all.
if [[ -n ${1:-} ]]; then
  CFG=$1
  [[ -f $CFG ]] || { echo "[failed] no such config: $CFG"; exit 1; }
else
  # With nothing named, look where a config actually ends up: here, then the
  # download folder. Reading each candidate is what makes this safe, so a JSON
  # that is not one of these is passed over rather than half-built.
  FOUND=$("$PY" "$ROOT/scripts/find_config.py" --explain 2>/dev/null)
  if [[ -z $FOUND ]]; then
    echo "[failed] no config named, and none found in $(pwd) or ~/Downloads."
    echo "         Run ./f127 new, or design one in docs/build.html and download it."
    exit 1
  fi
  CFG=${FOUND%%$'\t'*}
  echo "  reading $(basename "$CFG")  (${FOUND#*$'\t'})"
fi

read -r N_CHAINS BOX SALT SALT_M GUEST N_GUEST ROUTE SALTS N_INSIDE < <(
  "$PY" "$ROOT/scripts/read_config.py" "$CFG"
) || { echo "[failed] could not read $CFG"; exit 1; }
[[ $SALTS == "-" ]] && SALTS=""
export SALTS

[[ $N_CHAINS == 34 ]] || { echo "[failed] only 34 chains ship with the tool. Assembling another number is not implemented yet."; exit 1; }

# What this is going to cost, before it starts costing it. The rates are
# measured on the 577,146-atom loaded system, so they are the right order for
# anything this tool builds.
if [[ $RUN == 1 ]]; then
  GPU=no; GPU_NAME=""
  command -v gmx >/dev/null 2>&1 && gmx --version 2>/dev/null \
    | grep -qiE '^GPU support: *(CUDA|OpenCL|SYCL|HIP)' && GPU=yes
  [[ $GPU == yes ]] && GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
  "$PY" - "$NPT_NS" "$NS" "$GPU" "$GPU_NAME" <<'PLAN'
import sys
npt, prod, gpu, name = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3] == "yes", sys.argv[4].strip()
# ns/day on the 18 nm default box, as reported by people who ran the tool.
measured = [("5090", 106.3, "one RTX 5090"), ("3090", 43.1, "one RTX 3090"), ("3070", 39.4, "one RTX 3070")]
if gpu:
    rate, where = 43.1, "one RTX 3090"
    for key, r, w in measured:
        if key in name:
            rate, where = r, w
            break
    else:
        if name:
            where += f" (this machine has {name}, not measured yet)"
else:
    rate, where = 1.4, "10 CPU threads, no GPU"
h = lambda ns: ns / rate * 24
print(f"  plan: equilibrate {npt:g} ns, then produce {prod:g} ns")
print(f"        about {h(npt) + h(prod):.1f} h at {rate:g} ns/day, measured on {where}")
if npt > prod:
    print(f"        the equilibration is the longer of the two. --npt N changes it,")
    print(f"        and ./f127 test does a short pass of everything")
PLAN
fi

case $ROUTE in
  solution) METHODS="solution" ;;
  shell)    METHODS="shell" ;;
  both)     METHODS="solution shell" ;;
  *) echo "[failed] unknown route: $ROUTE"; exit 1 ;;
esac

for M in $METHODS; do
  OUT="run_${GUEST}_${M}"
  # the hollow of the shell template holds far fewer than the water around it,
  # so the two systems of a "both" run are counted separately
  N="$N_GUEST"
  [[ $M == shell ]] && N="${N_INSIDE:-$N_GUEST}"
  echo "=== $OUT  ($N x $GUEST, $M, box $BOX nm, ${SALTS:-$SALT $SALT_M M})"
  bash "$ROOT/scripts/2_load.sh" --guest "$GUEST" --n "$N" --method "$M" \
       --box "$BOX" --conc "$SALT_M" --out "$OUT" || exit 1
  [[ -s $OUT/ions.gro ]] || { echo "[failed] $OUT was not built"; exit 1; }
  cp "$CFG" "$OUT/$(basename "$CFG")"
  if [[ $RUN == 1 ]]; then
    # 2 ns of NpT is right on a GPU and most of a day on a laptop, so the
    # length is settable for anyone checking that the pipeline runs at all.
    bash "$ROOT/scripts/3_equilibrate.sh" "$OUT" "$NPT_NS" || exit 1
    bash "$ROOT/scripts/4_production.sh" "$OUT" "$NS" || exit 1
  else
    ( cd "$OUT" && ${GMX:-gmx} grompp -f "$ROOT/mdp/em.mdp" -c ions.gro -p topol.top \
        -n index.ndx -o em.tpr -maxwarn 5 >grompp_em.log 2>&1 ) \
      && echo "  built as far as em.tpr. Next: bash scripts/3_equilibrate.sh $OUT" \
      || { echo "  [failed] grompp. See $OUT/grompp_em.log"; exit 1; }
  fi
done
echo "done."
