#!/usr/bin/env bash
# bash scripts/4_production.sh run_dir 100      (ns)
set +u; source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true; set -uo pipefail
GMX=${GMX:-gmx}; ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; NS=${2:-100}; cd "$D"
sed "s/^nsteps.*/nsteps = $((NS*500000))/" "$ROOT/mdp/prod.mdp" > prod.mdp
$GMX grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o prod.tpr -maxwarn 5 >grompp_prod.log 2>&1 || exit 1
$GMX mdrun -deffnm prod ${MDRUN_OPT:--ntmpi 1 -ntomp 16 -nb gpu -bonded gpu -pme gpu -update gpu} >prod.out 2>&1
grep -i Performance prod.log | tail -1
