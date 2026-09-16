#!/usr/bin/env bash
# Make the archive that goes to another person.
#
#   bash scripts/9_ship.sh            -> ~/Downloads/f127-load_<date>.tar.gz
#   bash scripts/9_ship.sh /some/dir  -> writes it there instead
#
# Not the same thing as 9_pack.sh, which bundles one finished run. This bundles
# the tool.
#
# Left out: run_*/ (build output), compat.json (describes one machine),
# system.cfg and system.json (a config in the root would be picked up in place
# of a mistyped path), mdout.mdp, __pycache__, .git, .DS_Store and ._* files.
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
