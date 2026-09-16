#!/usr/bin/env bash
# Bundle a run so it can be moved, archived or handed over.
#
#   bash scripts/9_pack.sh run_pyrene_solution
#
# The bundle carries the config, the topology with its include files resolved,
# the coordinates, the index and the mdp files. It does not carry trajectories:
# those are large and are not what someone needs to reproduce the build.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}
[[ -d $D ]] || { echo "[failed] no $D"; exit 1; }
NAME=$(basename "$D")
STAGE=$(mktemp -d)
mkdir -p "$STAGE/$NAME"
# The repository's mdp files first, so em.mdp travels, then the run's own on top
# of them. The order used to be the other way round, which meant the bundle
# shipped the templates with nsteps unset and GUESTRESN unsubstituted rather than
# the files the run was actually built from.
cp "$ROOT"/mdp/*.mdp "$STAGE/$NAME/" 2>/dev/null
for f in system.cfg system.json topol.top index.ndx ions.gro em.tpr npt.gro prod.tpr \
         npt.mdp prod.mdp; do
  [[ -f $D/$f ]] && cp "$D/$f" "$STAGE/$NAME/"
done
# compat.json is written where ./f127 check was run, which is not the run
# directory, so it was listed above and never found.
[[ -f $ROOT/compat.json ]] && cp "$ROOT/compat.json" "$STAGE/$NAME/"
# toppar is a symlink into the repo and ff/ holds the guest topology together
# with the parameters CHARMM-GUI supplied for it. Both have to travel, or the
# bundle grompps on the machine that made it and nowhere else.
for sub in toppar ff; do
  [[ -d $D/$sub ]] || continue
  mkdir -p "$STAGE/$NAME/$sub"
  cp -L "$D/$sub"/*.itp "$STAGE/$NAME/$sub/" 2>/dev/null
done
cat > "$STAGE/$NAME/README.txt" <<TXT
Built with f127-load.  $(date -u '+%Y-%m-%d %H:%M UTC')

  system.cfg    the choices this system was built from
  topol.top     topology, include files are in toppar/
  toppar/       the polymer, water and ion parameters
  ff/           the solute topology and the parameters CHARMM-GUI supplied
  ions.gro      coordinates after solvation and ions
  index.ndx     System, micelle, solute, W_ION, core, corona, water, ions
  em.mdp        minimisation
  npt.mdp       equilibration, as this run used it
  prod.mdp      production, as this run used it

Trajectories are not in here. They are large and are not what is needed to
rebuild the system.

To continue:
  gmx grompp -f em.mdp -c ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5
  gmx mdrun -deffnm em
TXT
OUTF="$PWD/${NAME}.tar.gz"
tar -czf "$OUTF" -C "$STAGE" "$NAME"
rm -rf "$STAGE"
echo "  $OUTF  ($(du -h "$OUTF" | cut -f1))"
tar -tzf "$OUTF" | sed 's/^/    /'
