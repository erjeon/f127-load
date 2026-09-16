#!/usr/bin/env bash
# Build, equilibrate and run the two ends of what the tool accepts.
#
#   bash tests/edge_test.sh [ns]        default 3
#
# The least the tool will take is one solute molecule, no salt and a box far
# larger than the micelle. The most is the shell route filled to its limit with
# the four-salt physiological mixture.
#
# The small box is asked for but not obtained. The host micelle was equilibrated
# in 25 nm, so anything smaller is refused and built at 25 nm with a warning.
# Building a genuinely smaller box needs densify.py to compress the host first,
# which this test does not cover. Between them they exercise the placement, the ion mixing, the
# wrapping into a smaller box and the pressure coupling.
#
# This is the test that found three real defects: a failed stage that did not
# stop the run, a shell capacity that accepted more molecules than could be
# placed, and a GROMACS version being swapped underneath the caller.
set -uo pipefail
PY=${PYTHON:-python3}
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NS=${1:-3}
cd "$ROOT"

mk() {  # name, answers
  printf '%b' "$2" | "$PY" f127load/wizard.py >/dev/null 2>&1
  [[ -s $1 ]] || { echo "[failed] the wizard did not write $1"; exit 1; }
}
mk "$ROOT/edge_min.cfg" '34\n2\n30.0\n4\n5\n1\n1\n1\n'"$ROOT"'/edge_min.cfg\n'
# 28 rather than the 36 asked for before. The builder places at random and
# rejects overlaps, so it saturates at 32 for pyrene however many are asked for,
# and the wizard now stops below that. Asking for more only tests the refusal.
mk "$ROOT/edge_max.cfg" '34\n1\n14.5\n1\n3\n5\n1\n28\n2\n'"$ROOT"'/edge_max.cfg\n'

fail=0
for c in edge_min edge_max; do
  read -r _ _ _ _ GUEST _ ROUTE _ < <("$PY" scripts/read_config.py "$c.cfg")
  case $ROUTE in solution) M=solution;; shell) M=shell;; *) M="solution shell";; esac
  echo "########## $c  ($GUEST, $ROUTE)"
  ./f127 build "$c.cfg" || { echo "  [failed] build"; fail=1; continue; }
  for m in $M; do
    D="run_${GUEST}_${m}"
    bash scripts/3_equilibrate.sh "$D" 0.4 || { echo "  [failed] equilibrate"; fail=1; continue; }
    bash scripts/4_production.sh  "$D" "$NS" || { echo "  [failed] production"; fail=1; continue; }
    [[ -s $D/prod.gro ]] && echo "  [ok] $D ran $NS ns" || { echo "  [failed] $D"; fail=1; }
  done
done
echo "########## $([[ $fail == 0 ]] && echo "all passed" || echo "something failed")"
exit $fail
