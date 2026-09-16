#!/usr/bin/env bash
# bash scripts/4_production.sh run_dir 100      (ns)
source "$(dirname "${BASH_SOURCE[0]}")/_gmxenv.sh"
gmx_banner
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; NS=${2:-100}; cd "$D"
[[ -s npt.gro ]] || { echo "  [failed] no npt.gro in $D. Equilibrate first"; exit 1; }
# fractional ns are allowed, so the step count is computed in Python
STEPS=$("$PY" -c "import sys; print(round(float(sys.argv[1])*500000))" "$NS")
# The solute's residue name, read from the topology: anything that is not the
# polymer, water or an ion. The mdp files name a temperature-coupling group for it.
guest_resn() {
  awk '/^\[ *molecules/{f=1;next}
       f && NF==2 && $1!="S1P1" && $1!="TIP3" && $1!="SOL" && $1!="SOD" &&
       $1!="CLA" && $1!="POT" && $1!="CAL" && $1!="MG" {print $1; exit}' topol.top
}
RESN=$(guest_resn)
[[ -n ${RESN:-} ]] || { echo "  [failed] no solute found in topol.top"; exit 1; }

sed "s/^nsteps.*/nsteps = $STEPS/; s/GUESTRESN/$RESN/g${BAROSTAT_SED:+; }${BAROSTAT_SED}" "$ROOT/mdp/prod.mdp" > prod.mdp
echo "  $NS ns = $STEPS steps"
$GMX grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o prod.tpr -maxwarn 5 >grompp_prod.log 2>&1 \
  || { echo "  [failed] grompp before production. See $D/grompp_prod.log"
       tail -6 grompp_prod.log; exit 1; }
# -update gpu is left off, see 3_equilibrate.sh
$GMX mdrun -deffnm prod $(mdrun_opt md) >prod.out 2>&1
[[ -s prod.gro ]] || { echo "  [failed] production wrote no prod.gro. See prod.out"
                       tail -4 prod.out; exit 1; }
grep -i Performance prod.log | tail -1
