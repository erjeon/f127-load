#!/usr/bin/env bash
# Build a solute-loaded F127 micelle system.
#   bash scripts/2_load.sh --guest ibuprofen --n 8 --method shell --box 24.8 --out run_ibu
# Sourcing a GMXRC unconditionally put GROMACS 2022.3 in front of a 2025.4
# already on PATH, and trjconv then refused the 2025 tpr it was handed.
source "$(dirname "${BASH_SOURCE[0]}")/_gmxenv.sh"
gmx_line
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GUEST=""; N=8; METHOD=shell; BOX=""; OUT=""; CONC=0.154; SEED=2026
# Both forms work. The short one is what the documentation shows.
#   bash scripts/2_load.sh pyrene 100 solution
#   bash scripts/2_load.sh --guest pyrene --n 100 --method solution --out run_pyr
if [[ ${1:-} != --* && $# -ge 1 ]]; then
  GUEST=$1; [[ ${2:-} ]] && N=$2; [[ ${3:-} ]] && METHOD=$3; [[ ${4:-} ]] && OUT=$4
else
  while [[ $# -gt 0 ]]; do case $1 in
    --guest) GUEST=$2; shift 2;;  --n) N=$2; shift 2;;
    --method) METHOD=$2; shift 2;; --box) BOX=$2; shift 2;;
    --out) OUT=$2; shift 2;;      --conc) CONC=$2; shift 2;;
    --seed) SEED=$2; shift 2;;    *) echo "unknown: $1"; exit 1;; esac; done
fi
[[ -n $GUEST ]] || { echo "usage: 2_load.sh GUEST [N] [shell|solution] [OUTDIR]"; exit 1; }
# "solution" is what the paper and the documentation call it, "soak" is the old name
[[ $METHOD == solution ]] && METHOD=soak
OUT=${OUT:-run_${GUEST}_${METHOD}}
# the shell route closes a template that needs the larger box, the solution route
# starts from the equilibrated micelle and needs only its own
# Both routes start from a structure equilibrated in a 24.8 nm box, so that is
# the default for both. A smaller box is reached afterwards with densify.py,
# which removes water and lets the barostat compress, rather than by cutting the
# box around a structure that does not fit it.
# each host carries its own equilibrated box, and going below it folds the corona
[[ -n $BOX ]] || { [[ $METHOD == shell ]] && BOX=25.0 || BOX=24.8; }

LIB="$ROOT/library/$GUEST"
[[ -d $LIB ]] || { echo "[failed] no library/$GUEST"; exit 1; }
RESN=$(cat "$LIB/RESNAME.txt" 2>/dev/null | tr -d '[:space:]')
if [[ -z ${RESN:-} ]]; then
  # Six of the library folders hold a structure and nothing else. They are input
  # for parameter generation, not molecules that can be built, and saying only
  # that RESNAME.txt is missing sent people off to write one by hand.
  echo "[failed] $GUEST has no force field parameters yet."
  echo "  $LIB holds a structure to generate them from and nothing else."
  echo "  Put $LIB/$GUEST.mol2 through CHARMM-GUI Ligand Reader and Modeler,"
  echo "  drop what it returns into $LIB, then build the force field fragment:"
  echo "    python3 scripts/make_ff_fragment.py $LIB /path/to/charmm36.ff"
  echo "  library/README.md lists which molecules are ready."
  exit 1
fi
# ff_*.itp is the force field fragment, not the molecule. Picking the first
# .itp by name grabbed it whenever the residue name sorts after "ff", which
# is why paclitaxel broke and curcumin did not.
ITP=$(ls "$LIB"/*.itp 2>/dev/null | grep -v "/ff_" | head -1)
[[ -n ${ITP:-} ]] || { echo "[failed] no CGenFF itp in $LIB. Put the CHARMM-GUI output there."; exit 1; }
CRD=""
for c in "$LIB/${RESN}.pdb" "$LIB/${RESN}.gro" "$LIB/${GUEST}.pdb" "$LIB/${GUEST}.gro"; do
  [[ -f $c ]] && { CRD=$c; break; }
done
[[ -n $CRD ]] || { echo "[failed] no coordinates in $LIB. Put the ${RESN}.pdb CHARMM-GUI returned there."; exit 1; }
# The topology and the coordinates must describe the same atoms in the same order.
# A structure built from SMILES has the right count but generic names, so counting
# alone does not catch a mismatch. Compare the names.
"$PY" - "$ITP" "$CRD" <<'CHK' || exit 1
import re, sys
itp, crd = sys.argv[1], sys.argv[2]
sec = re.search(r"\[ *atoms *\](.*?)(?:\n\[|\Z)", open(itp).read(), re.S).group(1)
a = [l.split()[4] for l in sec.splitlines() if l.strip() and not l.lstrip().startswith(";")]
if crd.endswith(".pdb"):
    b = [l[12:16].strip() for l in open(crd) if l.startswith(("ATOM", "HETATM"))]
else:
    L = open(crd).read().splitlines(); b = [l[10:15].strip() for l in L[2:2 + int(L[1])]]
if len(a) != len(b):
    sys.exit(f"[failed] itp has {len(a)} atoms, {crd} has {len(b)}")
bad = [(i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
if bad:
    i, x, y = bad[0]
    sys.exit(f"[failed] atom {i+1} is {x} in the topology and {y} in {crd}. "
             "Use the coordinates CHARMM-GUI returned, not the structure you uploaded.")
print(f"  topology and coordinates agree on {len(a)} atoms")
CHK

# The residue name has to read the same in four places: the file that declares
# it, the moleculetype grompp looks up, the residue column of the topology and
# the residue column of the coordinates. Only the atom names were compared
# before, so a guest whose parameters came from another project under a
# different name built for five minutes and then stopped at grompp with
# "No such moleculetype". Doxorubicin went round eight times over exactly that.
"$PY" - "$RESN" "$ITP" "$CRD" <<'CHK' || exit 1
import re, sys
resn, itp, crd = sys.argv[1], sys.argv[2], sys.argv[3]
text = open(itp).read()

mt = re.search(r"\[ *moleculetype *\](.*?)(?:\n\[|\Z)", text, re.S).group(1)
declared = next(l.split()[0] for l in mt.splitlines()
                if l.strip() and not l.lstrip().startswith(";"))

sec = re.search(r"\[ *atoms *\](.*?)(?:\n\[|\Z)", text, re.S).group(1)
in_itp = {l.split()[3] for l in sec.splitlines()
          if l.strip() and not l.lstrip().startswith(";")}

if crd.endswith(".pdb"):
    in_crd = {l[17:20].strip() for l in open(crd) if l.startswith(("ATOM", "HETATM"))}
else:
    L = open(crd).read().splitlines()
    in_crd = {l[5:10].strip() for l in L[2:2 + int(L[1])]}

problems = []
if declared != resn:
    problems.append(f"RESNAME.txt says {resn}, the moleculetype is {declared}")
if in_itp != {resn}:
    problems.append(f"the topology labels its residues {sorted(in_itp)}")
if in_crd != {resn}:
    problems.append(f"{crd.rsplit('/', 1)[-1]} labels its residues {sorted(in_crd)}")
if problems:
    note = ""
    if len(resn) > 3:
        note = ("\n  A PDB residue name is three characters. A longer one runs into "
                "the chain column and is read back short.")
    sys.exit("[failed] the residue name has to match everywhere.\n  "
             + "\n  ".join(problems) + note)
print(f"  residue name {resn} agrees in the topology and the coordinates")
CHK

case $METHOD in
  # f127_shell_template.gro still carries the eight pyrene it was built with,
  # so loading a different guest into it would silently give a system holding
  # both. f127_shell_bare.gro is the same shell with those removed.
  shell) HOST="$ROOT/data/f127_shell_bare.gro";;
  soak)  HOST="$ROOT/data/f127_micelle_34.gro";;
  *) echo "method must be shell or solution"; exit 1;;
esac

mkdir -p "$OUT"; cd "$OUT"
ln -sfn "$ROOT/data/toppar" toppar; mkdir -p ff; cp "$ITP" ff/
# a guest usually needs atom types the shipped force field does not carry,
# because that file holds only what the polymer and water use
FF=$(ls "$LIB"/ff_*.itp 2>/dev/null | head -1)
[[ -n ${FF:-} ]] && cp "$FF" ff/
echo "[1/4] placing the solute ($METHOD, $N x $RESN)"
"$PY" "$ROOT/scripts/place_guest.py" "$HOST" "$CRD" "$N" "$METHOD" loaded.gro "$SEED" "$BOX" || exit 1
$GMX editconf -f loaded.gro -o boxed.gro -c -box $BOX $BOX $BOX >editconf.log 2>&1

# editconf sets the box and centres the contents, but it does not bring back
# anything that now lies outside. With a box smaller than the one the structure
# was equilibrated in, chain tips stick out by a couple of nm, solvate fills the
# space they will occupy once wrapped, and the first step of minimisation sees
# an infinite force. Folding them into the box first is what avoids that.
OUTSIDE=$("$PY" "$ROOT/scripts/wrap_into_box.py" boxed.gro boxed.gro) || exit 1
[[ $OUTSIDE == "0" ]] || echo "  wrapped $OUTSIDE atom(s) back into the box"

cat > topol.top <<TOP
#include "toppar/forcefield.itp"
$([[ -n ${FF:-} ]] && echo "#include \"ff/$(basename "$FF")\"")
#include "toppar/ion_types.itp"
#include "toppar/S1P1.itp"
#include "toppar/TIP3.itp"
#include "toppar/SOD.itp"
#include "toppar/CLA.itp"
#include "toppar/POT.itp"
#include "toppar/CAL.itp"
#include "toppar/MG.itp"
#include "ff/$(basename "$ITP")"
[ system ]
F127 micelle + $RESN
[ molecules ]
S1P1              34
$RESN             $N
TOP

echo "[2/4] solvating"
$GMX solvate -cp boxed.gro -cs "$ROOT/data/tip3p.gro" -o solv.gro -p topol.top >solvate.log 2>&1 \
  || { echo "[failed] solvate. See $(pwd)/solvate.log"; tail -6 solvate.log; exit 1; }
# solvate takes its radii from atom names, and a CHARMM-GUI name it does not
# recognise gets a guessed one that is too small. Paclitaxel came out with water
# 0.137 nm from a heavy atom and minimisation pulled it apart on the first step.
"$PY" "$ROOT/scripts/trim_close_water.py" solv.gro topol.top || exit 1
# Ions. SALTS is "name:cation:anion:n_cation:n_anion:cation_charge" entries
# separated by commas, written by the wizard. genion places one cation and one
# anion type per call and rewrites topol.top as it goes, so each salt needs its
# own grompp first.
#
# The charge has to be passed. genion defaults to -pq 1 and does not read the
# charge from the ion topology, so a divalent cation is counted as singly
# charged: asking for 8 Mg and 16 Cl gave 16 Mg and 16 Cl and a system at +16.
#
# Neutralising is a separate pass at the end with sodium and chloride, rather
# than -neutral hung on the last salt. That keeps the requested stoichiometry
# exact and puts any residual charge, which comes from the solute rather than
# from the salt, somewhere predictable.
if [[ -n ${SALTS:-} ]]; then
  IFS=',' read -ra ENTRIES <<< "$SALTS"
  IN=solv.gro
  for i in "${!ENTRIES[@]}"; do
    IFS=':' read -r SNAME SCAT SAN SNP SNN SPQ <<< "${ENTRIES[$i]}"
    SPQ=${SPQ:-1}
    echo "[3/4] ions, $SNAME ($SNP $SCAT at ${SPQ}+, $SNN $SAN)"
    $GMX grompp -f "$ROOT/mdp/em.mdp" -c "$IN" -p topol.top -o ions.tpr -maxwarn 5 \
      >"grompp_ions_$SNAME.log" 2>&1 || \
      { echo "[failed] grompp before $SNAME. See grompp_ions_$SNAME.log"; exit 1; }
    echo TIP3 | $GMX genion -s ions.tpr -o ions.gro -p topol.top \
      -pname "$SCAT" -nname "$SAN" -np "$SNP" -nn "$SNN" -pq "$SPQ" -nq -1 \
      >"genion_$SNAME.log" 2>&1 || \
      { echo "[failed] genion for $SNAME. See genion_$SNAME.log"; exit 1; }
    IN=ions.gro
  done
  echo "[3/4] ions, neutralising"
  $GMX grompp -f "$ROOT/mdp/em.mdp" -c "$IN" -p topol.top -o ions.tpr -maxwarn 5 \
    >grompp_ions_neutral.log 2>&1 || \
    { echo "[failed] grompp before neutralising. See grompp_ions_neutral.log"; exit 1; }
  echo TIP3 | $GMX genion -s ions.tpr -o ions.gro -p topol.top \
    -pname SOD -nname CLA -neutral >genion_neutral.log 2>&1 || \
    { echo "[failed] genion while neutralising. See genion_neutral.log"; exit 1; }
else
  echo "[3/4] ions ($CONC M NaCl)"
  $GMX grompp -f "$ROOT/mdp/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 5 >grompp_ions.log 2>&1 || \
    { echo "[failed] grompp for ions. See grompp_ions.log"; exit 1; }
  echo TIP3 | $GMX genion -s ions.tpr -o ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc $CONC >genion.log 2>&1 || exit 1
fi
echo "[4/4] index"
"$PY" "$ROOT/scripts/mkndx.py" ions.gro index.ndx || exit 1
tail -6 topol.top
# genion runs once per salt and GROMACS backs up ions.gro, ions.tpr and
# topol.top on every pass. A four-salt build left 130 MB of #file.N# behind.
# find -name '#*#' does not match these, the trailing # is part of the name but
# the shell still needs the leading one quoted.
rm -f ./\#*\# 2>/dev/null
echo "done: $(pwd)"
