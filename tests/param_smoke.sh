#!/usr/bin/env bash
# Check that a CGenFF parameter set actually runs.
#
# One molecule alone in a water box, minimised and then run briefly at constant
# pressure. This does not test whether the parameters are good, only that they
# are self-consistent enough for GROMACS to integrate. A set that fails here is
# broken. A set that passes may still be wrong, which is what check_params.py
# is for.
#
#   bash tests/param_smoke.sh <charmm-gui-dir> <RESNAME>
set +u
# a gmx already on PATH is the one the caller meant; see scripts/3_equilibrate.sh
command -v gmx >/dev/null 2>&1 || source "${GMXRC:-/usr/local/gromacs/bin/GMXRC}" 2>/dev/null || true
set -uo pipefail

GMX=${GMX:-gmx}
DIR=${1:?charmm-gui output directory}
RES=${2:?residue name, for example TAX}
NAME=$(basename "$DIR")
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT

fail() { printf "  failed  %s\n" "$*"; [ -f "$W/$3" ] && tail -12 "$W/$3"; exit 1; }

cp "$DIR/gromacs/$RES.itp" "$W/" 2>/dev/null || fail "no $RES.itp"
cp "$DIR/gromacs/charmm36.itp" "$W/" 2>/dev/null || fail "no charmm36.itp"
# charmm36.itp from CHARMM-GUI carries the force field but no water, so the
# TIP3 topology has to come from somewhere. Take the one shipped with the tool.
TOPPAR=${TOPPAR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/toppar}
cp "$TOPPAR/TIP3.itp" "$W/" 2>/dev/null || fail "no TIP3.itp in $TOPPAR"
# charmm36.itp as exported carries only the atom types the uploaded molecule
# uses, so the water types OT and HT are absent and have to be supplied
cp "$TOPPAR/water_types.itp" "$W/" 2>/dev/null || fail "no water_types.itp"
cp "$DIR/ligandrm.pdb" "$W/lig.pdb" 2>/dev/null || fail "no ligandrm.pdb"

cat > "$W/topol.top" <<TOP
#include "charmm36.itp"
#include "water_types.itp"
#include "$RES.itp"
#include "TIP3.itp"

[ system ]
$NAME in water

[ molecules ]
$RES   1
TOP

cd "$W"
printf "%-24s " "$NAME"

$GMX editconf -f lig.pdb -o box.gro -c -d 1.2 -bt cubic >e.log 2>&1 \
  || fail "editconf" "" e.log
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top >s.log 2>&1 \
  || fail "solvate" "" s.log
# gmx solvate writes SOL into the topology, but the CHARMM water moleculetype
# is called TIP3. Rename it in the topology only. The gro file is fixed format,
# so a plain substitution there shifts every column by one and grompp rejects
# the file.
sed -i.bak 's/^SOL/TIP3/' topol.top
NW=$(awk '$1=="TIP3"{print $2}' topol.top | tail -1)

cat > em.mdp <<'M'
integrator = steep
nsteps     = 2000
emtol      = 500
cutoff-scheme = Verlet
coulombtype   = PME
rvdw          = 1.2
rcoulomb      = 1.2
vdw-modifier  = force-switch
rvdw-switch   = 1.0
M
$GMX grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 3 >g1.log 2>&1 \
  || fail "grompp EM" "" g1.log
$GMX mdrun -deffnm em -ntmpi 1 -ntomp 4 -nb cpu >m1.log 2>&1 \
  || fail "EM" "" m1.log
FMAX=$(grep -m1 "Maximum force" em.log | awk '{print $4}')

cat > npt.mdp <<'M'
integrator  = md
dt          = 0.002
nsteps      = 25000
nstlog      = 5000
nstenergy   = 500
nstxout     = 0
nstvout     = 0
nstfout     = 0
cutoff-scheme = Verlet
coulombtype   = PME
rvdw          = 1.2
rcoulomb      = 1.2
vdw-modifier  = force-switch
rvdw-switch   = 1.0
constraints   = h-bonds
tcoupl        = v-rescale
tc_grps       = System
tau_t         = 1.0
ref_t         = 310.15
pcoupl        = C-rescale
tau_p         = 5.0
compressibility = 4.5e-5
ref_p         = 1.0
gen-vel       = yes
gen-temp      = 310.15
M
$GMX grompp -f npt.mdp -c em.gro -p topol.top -o npt.tpr -maxwarn 3 >g2.log 2>&1 \
  || fail "grompp NPT" "" g2.log
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp 4 -nb cpu >m2.log 2>&1 \
  || fail "NPT 50 ps" "" m2.log

# gmx energy prints its summary table to stderr, so read the trajectory of the
# temperature out of the xvg it writes and average the second half
printf "Temperature\n" | $GMX energy -f npt.edr -o temp.xvg >/dev/null 2>&1
T=$(awk '!/^[#@]/{v[n++]=$2} END{s=0; for(i=int(n/2);i<n;i++) s+=v[i];
         if(n>1) printf "%.1f", s/(n-int(n/2)); else print 0}' temp.xvg)
LINCS=$(grep -c "LINCS WARNING" npt.log 2>/dev/null | head -1)
CONSTR=$(awk '/Constr. rmsd/{getline; print $NF}' npt.log | tail -1)

VERDICT=ok
awk -v t="${T:-0}" 'BEGIN{exit !(t>295 && t<325)}' || VERDICT="temperature out of range"
[ "${LINCS:-0}" -gt 0 ] 2>/dev/null && VERDICT="LINCS warnings"
printf "atoms %5s  water %5s  Fmax %11s  T %6.1f K  LINCS %2s  rmsd %10s  %s\n" \
  "$(sed -n 2p solv.gro | tr -d ' ')" "$NW" "$FMAX" "${T:-0}" "${LINCS:-0}" \
  "${CONSTR:-?}" "$VERDICT"
