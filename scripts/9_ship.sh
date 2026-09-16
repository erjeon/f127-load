#!/usr/bin/env bash
# Make the archive that goes to another person.
#
#   bash scripts/9_ship.sh            -> ~/Downloads/f127-load_<date>.tar.gz
#   bash scripts/9_ship.sh /some/dir  -> writes it there instead
#
# Not the same thing as 9_pack.sh, which bundles one finished run. This bundles
# the tool.
#
# What is left out and why:
#   run_*/          a build's output. One of them is 300 MB.
#   compat.json     written by ./f127 check, and it describes this machine
#   system.cfg
#   system.json     ★ a config left in the root is picked up when the path given
#                   on the command line does not exist, so shipping one means
#                   the person builds a system they did not ask for. That is
#                   exactly what happened here.
#   mdout.mdp       grompp leaves it behind
#   __pycache__     and .pyc
#   .git            the history is on GitHub, and it is larger than the tool
#   .DS_Store ._*   macOS
#
# Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR=${1:-$HOME/Downloads}
NAME="f127-load_$(date +%y%m%d)"
OUT="$OUTDIR/$NAME.tar.gz"

cd "$ROOT/.."
tar --exclude="$(basename "$ROOT")/run_*" \
    --exclude="$(basename "$ROOT")/.git" \
    --exclude="*/__pycache__" --exclude="*.pyc" \
    --exclude="*/compat.json" --exclude="*/system.cfg" \
    --exclude="*/system.json" --exclude="*/mdout.mdp" \
    --exclude=".DS_Store" --exclude="._*" \
    --exclude="*.tar.gz" \
    -czf "$OUT" --disable-copyfile "$(basename "$ROOT")" 2>/dev/null \
  || tar --exclude="$(basename "$ROOT")/run_*" \
         --exclude="$(basename "$ROOT")/.git" \
         --exclude="*/__pycache__" --exclude="*.pyc" \
         --exclude="*/compat.json" --exclude="*/system.cfg" \
         --exclude="*/system.json" --exclude="*/mdout.mdp" \
         --exclude=".DS_Store" --exclude="._*" --exclude="*.tar.gz" \
         -czf "$OUT" "$(basename "$ROOT")"

echo "  $OUT  ($(du -h "$OUT" | cut -f1))"
echo "  $(tar -tzf "$OUT" | wc -l | tr -d ' ') files"
# anything that should not have travelled
BAD=$(tar -tzf "$OUT" | grep -E "run_|compat\.json|system\.(cfg|json)|mdout\.mdp|__pycache__|\.git/" || true)
if [[ -n $BAD ]]; then
  echo "  [failed] these should not be in the archive:"; printf '%s\n' "$BAD" | sed 's/^/      /'; exit 1
fi
echo "  nothing that should have been left out is in it"
