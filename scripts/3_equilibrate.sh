#!/usr/bin/env bash
# EM → NPT.  bash scripts/3_equilibrate.sh run_dir [ns]
set +u; source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true; set -uo pipefail
GMX=${GMX:-gmx}; ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; NS=${2:-2}; cd "$D"
$GMX grompp -f "$ROOT/mdp/em.mdp" -c ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5 >grompp_em.log 2>&1 || exit 1
$GMX mdrun -deffnm em ${MDRUN_OPT:--ntmpi 1 -ntomp 16 -nb gpu} >em.out 2>&1
grep -m1 "Potential Energy" em.log; grep -m1 "Maximum force" em.log
FMAX=$(grep -m1 "Maximum force" em.log | awk '{print $4}')
awk -v f="$FMAX" 'BEGIN{ if (f+0 > 5000) print "  [경고] Fmax 가 크다. 겹침이 남아 있을 수 있다." }'
sed "s/SEEDVAL/$RANDOM/; s/^nsteps.*/nsteps = $((NS*500000))/" "$ROOT/mdp/npt.mdp" > npt.mdp
$GMX grompp -f npt.mdp -c em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 5 >grompp_npt.log 2>&1 || exit 1
$GMX mdrun -deffnm npt ${MDRUN_OPT:--ntmpi 1 -ntomp 16 -nb gpu -bonded gpu -pme gpu -update gpu} >npt.out 2>&1
[[ -s npt.gro ]] && { echo -n "  최종 box: "; tail -1 npt.gro; } || echo "  [FAIL] npt"
