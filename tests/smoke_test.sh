#!/usr/bin/env bash
# Checks only that the parameters run under GROMACS. It says nothing about loading.
#   bash tests/smoke_test.sh ibuprofen
# Use the gmx that is already on the PATH. Anyone running MD has set one up, and
# choosing a different one for them is not this script's business: a version is
# sometimes pinned on purpose. A GMXRC is sourced only when there is no gmx at
# all. Which one is in use gets printed, because the failure that started this
# was a silent swap from 2025.4 to 2022.3.
set +u
command -v gmx >/dev/null 2>&1 || source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true
set -uo pipefail
GMX=${GMX:-gmx}
echo "  using $($GMX --version 2>/dev/null | grep -m1 -o 'GROMACS.*' || echo "gmx (version unknown)")"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
G=${1:?guest name}; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
echo "=== smoke test: $G ==="
bash "$ROOT/scripts/2_load.sh" --guest "$G" --n 2 --method shell --box 22 --out "$T/s" >"$T/load.log" 2>&1 \
  || { echo "  [failed] building the system"; tail -15 "$T/load.log"; exit 1; }
echo "  [ok] built"
cd "$T/s"
sed 's/^nsteps.*/nsteps = 2000/' "$ROOT/mdp/em.mdp" > em_s.mdp
$GMX grompp -f em_s.mdp -c ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5 >g1.log 2>&1 \
  || { echo "  [failed] grompp"; tail -15 g1.log; exit 1; }
echo "  [ok] grompp"
$GMX mdrun -deffnm em ${MDRUN_OPT:--ntmpi 1 -ntomp 8 -nb cpu} >m1.log 2>&1 || { echo "  [failed] EM"; tail -15 m1.log; exit 1; }
FM=$(grep -m1 "Maximum force" em.log | awk '{print $4}')
echo "  [ok] EM  Fmax = $FM"
cat > md_s.mdp <<'M'
integrator = md
dt = 0.001
nsteps = 2000
cutoff-scheme = Verlet
nstlist = 20
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw_switch = 1.0
rvdw = 1.2
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
tcoupl = v-rescale
tc_grps = System
tau_t = 1.0
ref_t = 310.15
constraints = h-bonds
gen-vel = yes
gen-temp = 310.15
M
$GMX grompp -f md_s.mdp -c em.gro -p topol.top -n index.ndx -o md.tpr -maxwarn 5 >g2.log 2>&1 \
  || { echo "  [failed] grompp(md)"; tail -15 g2.log; exit 1; }
$GMX mdrun -deffnm md ${MDRUN_OPT:--ntmpi 1 -ntomp 8 -nb cpu} >m2.log 2>&1 \
  || { echo "  [failed] 2 ps of MD will not run"; tail -20 md.log; exit 1; }
echo "  [ok] 2 ps of MD completed"
echo "=== passed: $G ==="
