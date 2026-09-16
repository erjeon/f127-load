#!/usr/bin/env bash
# The eight standard analyses, plus figures.  bash scripts/5_analyze.sh run_dir [--begin ps]
# --begin must sit after equilibration. Averaging the whole trajectory mixes the
# collapsing structure with the equilibrated one and puts a hole in the core density.
source "$(dirname "${BASH_SOURCE[0]}")/_gmxenv.sh"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; shift; BEG=0
while [[ $# -gt 0 ]]; do case $1 in --begin) BEG=$2; shift 2;; *) shift;; esac; done
cd "$D"; O=results; mkdir -p $O
gmx_line
# check the imports before the pbc step, which takes minutes
"$PY" -c "import MDAnalysis, matplotlib" 2>/dev/null || {
  echo "  [failed] $PY cannot import MDAnalysis and matplotlib."
  echo "           The simulation is fine, only the analysis needs them."
  echo "           pip install \"MDAnalysis>=2.8\" matplotlib   in a virtual environment,"
  echo "           or point at one you have:  PYTHON=/path/to/venv/bin/python"
  exit 1
}
# the solute is whatever is not the polymer, water or an ion
RESN=$(awk '/^\[ *molecules/{f=1;next}
            f && NF==2 && $1!="S1P1" && $1!="TIP3" && $1!="SOL" && $1!="SOD" &&
            $1!="CLA" && $1!="POT" && $1!="CAL" && $1!="MG" {print $1; exit}' topol.top)
echo "  solute = $RESN, begin = $BEG ps"

if [[ ! -s $O/proc.xtc ]]; then
  echo "[1/8] periodic boundary"
  # everything downstream reads proc.xtc, so a failure here stops the run
  echo System | $GMX trjconv -s prod.tpr -f prod.xtc -n index.ndx -pbc whole -o $O/whole.xtc >$O/pbc1.log 2>&1 \
    || { echo "  [failed] trjconv -pbc whole. See $O/pbc1.log"; exit 1; }
  printf 'micelle\nSystem\n' | $GMX trjconv -s prod.tpr -f $O/whole.xtc -n index.ndx -pbc mol -center -o $O/proc.xtc >$O/pbc2.log 2>&1 \
    || { echo "  [failed] trjconv -pbc mol -center. See $O/pbc2.log"; exit 1; }
  rm -f $O/whole.xtc
fi
[[ -s $O/proc.xtc ]] || { echo "  [failed] $O/proc.xtc was not written"; exit 1; }
# core, corona and water are groups written by scripts/mkndx.py. An index given
# with -n replaces the moleculetype names, so TIP3 cannot be used here.
echo "[2/8] Rg"
echo micelle | $GMX gyrate -s prod.tpr -f $O/proc.xtc -n index.ndx -b $BEG -o $O/rg.xvg >$O/rg.log 2>&1 \
  || echo "  [failed] gyrate. See $D/$O/rg.log"
echo "[3/8] SASA"
echo micelle | $GMX sasa -s prod.tpr -f $O/proc.xtc -n index.ndx -b $BEG -o $O/sasa.xvg -dt 200 >$O/sasa.log 2>&1 \
  || echo "  [failed] sasa. See $D/$O/sasa.log"
echo "[4/8] RDF"
$GMX rdf -s prod.tpr -f $O/proc.xtc -n index.ndx -ref $RESN -sel core corona water -b $BEG -bin 0.02 -o $O/rdf.xvg >$O/rdf.log 2>&1 \
  || { echo "  [failed] rdf. See $D/$O/rdf.log"
       grep -A3 -m1 -iE "fatal|inconsistency" $O/rdf.log | sed 's/^/      /'; }
echo "[5/8] radial density, uptake, hydration"
"$PY" "$ROOT/scripts/radial.py" prod.tpr $O/proc.xtc $RESN $BEG $O || echo "  (radial failed)"
echo "[6/8] solute-solute association"
# one solute cannot associate with anything, which is a fact about the system
# and not a failure of the analysis
"$PY" "$ROOT/scripts/pipi.py" prod.tpr $O/proc.xtc $RESN $BEG $O \
  || echo "  (skipped, see the line above)"
echo "[7/8] interaction energy, by rerun"
if [[ ${DO_LIE:-1} == 1 ]]; then
  cat > $O/rerun.mdp <<'M'
integrator = md
nsteps = 0
cutoff-scheme = Verlet
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw_switch = 1.0
rvdw = 1.2
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
constraints = h-bonds
M
  echo "energygrps = $RESN core water" >> $O/rerun.mdp
  $GMX grompp -f $O/rerun.mdp -c npt.gro -p topol.top -n index.ndx -o $O/rerun.tpr -maxwarn 5 >$O/grompp_rerun.log 2>&1 \
   && $GMX mdrun -s $O/rerun.tpr -rerun $O/proc.xtc -e $O/rerun.edr -g $O/rerun.log -nb cpu -ntmpi 1 -ntomp $NTOMP >$O/rerun.out 2>&1 \
   && printf "Coul-SR:$RESN-core\nLJ-SR:$RESN-core\nCoul-SR:$RESN-water\nLJ-SR:$RESN-water\n\n" \
      | $GMX energy -f $O/rerun.edr -o $O/lie.xvg >$O/energy.log 2>&1 || echo "  (LIE failed)"
fi
echo "[8/8] figures"
"$PY" "$ROOT/scripts/plot_all.py" $O || echo "  (plot failed)"
ls -la $O/*.xvg $O/*.dat $O/figures/*.png 2>/dev/null | tail -20
