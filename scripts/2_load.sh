#!/usr/bin/env bash
# 화물이 담지된 F127 미셀 계를 만든다.
#   bash scripts/2_load.sh --guest ibuprofen --n 8 --method shell --box 24.8 --out run_ibu
set +u; source ${GMXRC:-/usr/local/gromacs/bin/GMXRC} 2>/dev/null || true; set -uo pipefail
GMX=${GMX:-gmx}
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GUEST=""; N=8; METHOD=shell; BOX=24.8; OUT=""; CONC=0.154; SEED=2026
while [[ $# -gt 0 ]]; do case $1 in
  --guest) GUEST=$2; shift 2;;  --n) N=$2; shift 2;;
  --method) METHOD=$2; shift 2;; --box) BOX=$2; shift 2;;
  --out) OUT=$2; shift 2;;      --conc) CONC=$2; shift 2;;
  --seed) SEED=$2; shift 2;;    *) echo "unknown: $1"; exit 1;; esac; done
[[ -n $GUEST && -n $OUT ]] || { echo "usage: --guest NAME --out DIR [--n 8] [--method shell|soak] [--box 24.8]"; exit 1; }

LIB="$ROOT/library/$GUEST"
[[ -d $LIB ]] || { echo "[FAIL] library/$GUEST 없음"; exit 1; }
RESN=$(cat "$LIB/RESNAME.txt" 2>/dev/null | tr -d '[:space:]')
[[ -n ${RESN:-} ]] || { echo "[FAIL] $LIB/RESNAME.txt 가 없다. library/README.md 참고."; exit 1; }
ITP=$(ls "$LIB"/*.itp 2>/dev/null | head -1)
[[ -n ${ITP:-} ]] || { echo "[FAIL] $LIB 에 CGenFF itp 가 없다. CHARMM-GUI 출력을 넣어라."; exit 1; }
CRD="$LIB/${RESN}.pdb"; [[ -f $CRD ]] || CRD="$LIB/${RESN}.gro"
[[ -f $CRD ]] || { echo "[FAIL] CHARMM-GUI 좌표 ${RESN}.pdb 가 없다 (itp 와 원자 순서가 같아야 한다)."; exit 1; }

case $METHOD in
  shell) HOST="$ROOT/data/f127_shell_template.gro";;
  soak)  HOST="$ROOT/data/f127_micelle_34.gro";;
  *) echo "method 는 shell 또는 soak"; exit 1;;
esac

mkdir -p "$OUT"; cd "$OUT"
ln -sfn "$ROOT/data/toppar" toppar; cp "$ITP" toppar/
echo "[1/4] 화물 배치 ($METHOD, $N x $RESN)"
python3 "$ROOT/scripts/place_guest.py" "$HOST" "$CRD" "$N" "$METHOD" loaded.gro "$SEED" || exit 1
$GMX editconf -f loaded.gro -o boxed.gro -c -box $BOX $BOX $BOX >editconf.log 2>&1

cat > topol.top <<TOP
#include "toppar/forcefield.itp"
#include "toppar/S1P1.itp"
#include "toppar/TIP3.itp"
#include "toppar/SOD.itp"
#include "toppar/CLA.itp"
#include "toppar/$(basename "$ITP")"
[ system ]
F127 micelle + $RESN
[ molecules ]
S1P1              34
$RESN             $N
TOP

echo "[2/4] 용매화"
$GMX solvate -cp boxed.gro -cs "$ROOT/data/tip3p.gro" -o solv.gro -p topol.top >solvate.log 2>&1 || exit 1
echo "[3/4] 이온 ($CONC M NaCl)"
$GMX grompp -f "$ROOT/mdp/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 5 >grompp_ions.log 2>&1 || \
  { echo "[FAIL] grompp(ions) — grompp_ions.log 확인"; exit 1; }
echo TIP3 | $GMX genion -s ions.tpr -o ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc $CONC >genion.log 2>&1 || exit 1
echo "[4/4] 인덱스"
python3 "$ROOT/scripts/mkndx.py" ions.gro index.ndx || exit 1
tail -6 topol.top
echo "완료: $(pwd)"
