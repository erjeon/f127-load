#!/usr/bin/env bash
# 파라미터가 GROMACS 에서 실제로 도는지만 확인한다. 담지 여부는 보지 않는다.
#   bash tests/smoke_test.sh ibuprofen
set +u; source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true; set -uo pipefail
GMX=${GMX:-gmx}; ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
G=${1:?guest name}; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
echo "=== smoke test: $G ==="
bash "$ROOT/scripts/2_load.sh" --guest "$G" --n 2 --method shell --box 22 --out "$T/s" >"$T/load.log" 2>&1 \
  || { echo "  [FAIL] 시스템 빌드"; tail -15 "$T/load.log"; exit 1; }
echo "  [ok] 빌드"
cd "$T/s"
sed 's/^nsteps.*/nsteps = 2000/' "$ROOT/mdp/em.mdp" > em_s.mdp
$GMX grompp -f em_s.mdp -c ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5 >g1.log 2>&1 \
  || { echo "  [FAIL] grompp"; tail -15 g1.log; exit 1; }
echo "  [ok] grompp"
$GMX mdrun -deffnm em ${MDRUN_OPT:--ntmpi 1 -ntomp 8 -nb cpu} >m1.log 2>&1 || { echo "  [FAIL] EM"; tail -15 m1.log; exit 1; }
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
  || { echo "  [FAIL] grompp(md)"; tail -15 g2.log; exit 1; }
$GMX mdrun -deffnm md ${MDRUN_OPT:--ntmpi 1 -ntomp 8 -nb cpu} >m2.log 2>&1 \
  || { echo "  [FAIL] 2 ps MD 가 돌지 않는다"; tail -20 md.log; exit 1; }
echo "  [ok] 2 ps MD 완주"
echo "=== PASS: $G ==="
