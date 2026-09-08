#!/usr/bin/env bash
# 표준 분석 8종 + PNG.   bash scripts/5_analyze.sh run_dir [--begin ps]
# ★ --begin 은 반드시 평형 도달 이후로 줄 것. 전 구간 평균은 코어 밀도에 인공물을 만든다.
set +u; source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true; set -uo pipefail
GMX=${GMX:-gmx}; ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=${1:?run_dir}; shift; BEG=0
while [[ $# -gt 0 ]]; do case $1 in --begin) BEG=$2; shift 2;; *) shift;; esac; done
cd "$D"; O=results; mkdir -p $O
RESN=$(awk '/^\[ molecules/{f=1;next} f&&NF==2&&$1!="S1P1"&&$1!="TIP3"&&$1!="SOD"&&$1!="CLA"{print $1;exit}' topol.top)
echo "  guest = $RESN, begin = $BEG ps"

if [[ ! -s $O/proc.xtc ]]; then
  echo "[1/8] PBC 보정"
  echo System | $GMX trjconv -s prod.tpr -f prod.xtc -n index.ndx -pbc whole -o $O/whole.xtc >$O/pbc1.log 2>&1
  printf 'micelle\nSystem\n' | $GMX trjconv -s prod.tpr -f $O/whole.xtc -n index.ndx -pbc mol -center -o $O/proc.xtc >$O/pbc2.log 2>&1
  rm -f $O/whole.xtc
fi
echo "[2/8] Rg";   echo micelle | $GMX gyrate -s prod.tpr -f $O/proc.xtc -n index.ndx -b $BEG -o $O/rg.xvg   >$O/rg.log 2>&1
echo "[3/8] SASA"; echo micelle | $GMX sasa   -s prod.tpr -f $O/proc.xtc -n index.ndx -b $BEG -o $O/sasa.xvg -dt 200 >$O/sasa.log 2>&1
echo "[4/8] RDF";  $GMX rdf -s prod.tpr -f $O/proc.xtc -n index.ndx -ref $RESN -sel core corona TIP3 -b $BEG -bin 0.02 -o $O/rdf.xvg >$O/rdf.log 2>&1
echo "[5/8] 반경 밀도 · 담지율 · 수화수"
python3 "$ROOT/scripts/radial.py" prod.tpr $O/proc.xtc $RESN $BEG $O || echo "  (radial 실패)"
echo "[6/8] guest-guest 회합"
python3 "$ROOT/scripts/pipi.py" prod.tpr $O/proc.xtc $RESN $BEG $O || echo "  (pipi 실패)"
echo "[7/8] 상호작용 에너지 (rerun)"
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
  echo "energygrps = $RESN core TIP3" >> $O/rerun.mdp
  $GMX grompp -f $O/rerun.mdp -c npt.gro -p topol.top -n index.ndx -o $O/rerun.tpr -maxwarn 5 >$O/grompp_rerun.log 2>&1 \
   && $GMX mdrun -s $O/rerun.tpr -rerun $O/proc.xtc -e $O/rerun.edr -g $O/rerun.log -nb cpu -ntmpi 1 -ntomp 16 >$O/rerun.out 2>&1 \
   && printf "Coul-SR:$RESN-core\nLJ-SR:$RESN-core\nCoul-SR:$RESN-TIP3\nLJ-SR:$RESN-TIP3\n\n" \
      | $GMX energy -f $O/rerun.edr -o $O/lie.xvg >$O/energy.log 2>&1 || echo "  (LIE 실패)"
fi
echo "[8/8] 그림"
python3 "$ROOT/scripts/plot_all.py" $O || echo "  (plot 실패)"
ls -la $O/*.xvg $O/*.dat $O/figures/*.png 2>/dev/null | tail -20
