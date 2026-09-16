#!/usr/bin/env bash
# EM → NPT.  bash scripts/3_equilibrate.sh run_dir [ns]
source "$(dirname "${BASH_SOURCE[0]}")/_gmxenv.sh"
gmx_banner
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; NS=${2:-2}; cd "$D"
$GMX grompp -f "$ROOT/mdp/em.mdp" -c ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5 >grompp_em.log 2>&1 \
  || { echo "  [failed] grompp before minimisation. See $D/grompp_em.log"
       tail -6 grompp_em.log; exit 1; }
# mdrun failing was not stopping anything. Minimisation would blow up, the next
# stage would grompp against a file that was never written, and the run would go
# on to production and report success. Every stage now checks its own output.
$GMX mdrun -deffnm em $(mdrun_opt em) >em.out 2>&1
[[ -s em.gro ]] || { echo "  [failed] minimisation wrote no em.gro. See em.out"
                     tail -4 em.out; exit 1; }
# mdrun dumps the whole system to step<N>b.pdb and step<N>c.pdb whenever LINCS
# complains, which it does for the first few dozen steps of a freshly packed box.
# Six of those left 262 MB of diagnostics beside a minimisation that converged.
# They are only worth keeping when it did not, so they go after em.gro exists.
rm -f step*b.pdb step*c.pdb em.trr
grep -m1 "Potential Energy" em.log; grep -m1 "Maximum force" em.log
FMAX=$(grep -m1 "Maximum force" em.log | awk '{print $4}')
awk -v f="$FMAX" 'BEGIN{ if (f+0 > 5000) print "  [warn] Fmax is large. Something may still overlap." }'
# The solute's residue name, read from the topology the build wrote. The mdp
# files name a temperature-coupling group for it, and that name used to be PYR,
# so every guest except pyrene stopped at grompp. Anything that is not the
# polymer, water or an ion is the solute.
guest_resn() {
  awk '/^\[ *molecules/{f=1;next}
       f && NF==2 && $1!="S1P1" && $1!="TIP3" && $1!="SOL" && $1!="SOD" &&
       $1!="CLA" && $1!="POT" && $1!="CAL" && $1!="MG" {print $1; exit}' topol.top
}
RESN=$(guest_resn)
[[ -n ${RESN:-} ]] || { echo "  [failed] no solute found in topol.top"; exit 1; }
echo "  solute $RESN"

# fractional ns, so a quick check does not have to run the full two
STEPS=$("$PY" -c "import sys; print(round(float(sys.argv[1])*500000))" "$NS")
sed "s/SEEDVAL/$RANDOM/; s/^nsteps.*/nsteps = $STEPS/; s/GUESTRESN/$RESN/g${BAROSTAT_SED:+; }${BAROSTAT_SED}" \
    "$ROOT/mdp/npt.mdp" > npt.mdp
echo "  npt $NS ns = $STEPS steps"
$GMX grompp -f npt.mdp -c em.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 5 >grompp_npt.log 2>&1 \
  || { echo "  [failed] grompp before NpT. See $D/grompp_npt.log"
       tail -6 grompp_npt.log; exit 1; }
# -update gpu is left off. It gives CUDA error #700, an illegal memory access,
# on this system during the first pressure-coupled steps. Constraints and the
# update run on the CPU instead, which costs a little and finishes.
$GMX mdrun -deffnm npt $(mdrun_opt md) >npt.out 2>&1
[[ -s npt.gro ]] || { echo "  [failed] npt wrote no npt.gro. See npt.out"
                      tail -4 npt.out; exit 1; }
echo -n "  final box: "; tail -1 npt.gro
