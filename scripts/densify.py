#!/usr/bin/env python3
"""Remove water and ions at random to reach a target polymer weight fraction.

The box is left alone and the barostat does the compressing. Shrinking the box
around a fixed configuration overlaps chains instead.

    python densify.py in.gro out.gro topol.top 0.10 [salt_M]

Everything that is not water or an ion is kept, so the micelle and whatever it
carries survive whatever their atom count. An earlier version took the first
70,010 atoms as the micelle, which was the polymer plus exactly eight pyrene.
With any other number of solute molecules it cut the rest away without a word:
28 solutes lost 20 of them.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import random
import sys
from pathlib import Path

SOLVENT = {"TIP3", "SOL", "HOH", "WAT"}
CATION = {"SOD", "NA", "POT", "K", "CAL", "CA", "MG", "ZN"}
ANION = {"CLA", "CL", "BR", "IOD"}
IONS = CATION | ANION

MW_F127 = 12586.0
N_CHAINS = 34
MW_W = 18.015
NA = 6.02214076e23


def usage():
    sys.exit(__doc__)


def residues(atoms):
    """Group .gro lines into residues, keyed by the residue number field."""
    out, cur, prev = [], [], None
    for a in atoms:
        rid = a[:5]
        if rid != prev and cur:
            out.append(cur)
            cur = []
        prev = rid
        cur.append(a)
    if cur:
        out.append(cur)
    return out


def resname(res):
    return res[0][5:10].strip()


def main():
    if len(sys.argv) not in (5, 6):
        usage()
    src, dst, top = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    wt = float(sys.argv[4])
    salt_m = float(sys.argv[5]) if len(sys.argv) == 6 else 0.154
    if not 0.0 < wt < 1.0:
        sys.exit(f"  [failed] weight fraction must be between 0 and 1, got {wt}")
    random.seed(20260903)

    lines = src.read_text().splitlines()
    n = int(lines[1])
    head, atoms, box = lines[0], lines[2:2 + n], lines[2 + n]
    if len(atoms) != n:
        sys.exit(f"  [failed] {src} says {n} atoms and carries {len(atoms)}")

    res = residues(atoms)
    keep_always = [r for r in res if resname(r) not in SOLVENT | IONS]
    wat = [r for r in res if resname(r) in SOLVENT]
    cats = [r for r in res if resname(r) in CATION]
    ans = [r for r in res if resname(r) in ANION]
    solute = sorted({resname(r) for r in keep_always})
    n_kept = sum(len(r) for r in keep_always)
    print(f"  keeping  {n_kept} atoms in {len(keep_always)} residues: "
          f"{', '.join(solute)}")
    print(f"  before   water {len(wat)}  cations {len(cats)}  anions {len(ans)}")

    m_pol = N_CHAINS * MW_F127 / NA
    n_wat = int(round(m_pol * (1 - wt) / wt / (MW_W / NA)))
    v_l = m_pol / wt * 1e-3
    n_ion = int(round(salt_m * v_l * NA))
    if n_wat > len(wat):
        # water is only ever removed, so a target below the box's present
        # weight fraction cannot be reached. The message used to say the
        # opposite, that the target was too dense.
        here = N_CHAINS * MW_F127 / (N_CHAINS * MW_F127 + len(wat) * MW_W)
        sys.exit(f"  [failed] {wt*100:.2f} wt% needs {n_wat} water molecules "
                 f"and the box holds {len(wat)}. This box is already at "
                 f"{here*100:.2f} wt%, and densify only removes water, so it "
                 f"cannot be made more dilute. Build in a larger box instead.")
    print(f"  target   {wt*100:.1f} wt%, {salt_m} M salt -> water {n_wat}, "
          f"{n_ion} of each ion")
    # one litre is 1e24 cubic nanometres. The figure printed here was three
    # orders out, which made a 19 nm box read as 1.9 nm.
    print(f"  the barostat should settle near "
          f"{(v_l * 1e24) ** (1 / 3):.2f} nm")

    out = list(keep_always)
    out += [wat[i] for i in sorted(random.sample(range(len(wat)), n_wat))]
    kc = sorted(random.sample(range(len(cats)), min(n_ion, len(cats))))
    ka = sorted(random.sample(range(len(ans)), min(n_ion, len(ans))))
    out += [cats[i] for i in kc]
    out += [ans[i] for i in ka]

    written, rid, aid = [], 0, 0
    for r in out:
        rid += 1
        for a in r:
            aid += 1
            written.append(f"{rid % 100000:5d}{a[5:15]}{aid % 100000:5d}{a[20:]}")
    dst.write_text(f"{head}\n{len(written)}\n" + "\n".join(written) + f"\n{box}\n")

    # the counts in the topology have to follow, one line per molecule name
    counts = {}
    for r in out:
        counts[resname(r)] = counts.get(resname(r), 0) + 1
    new = []
    for line in top.read_text().splitlines():
        k = line.split()
        if len(k) == 2 and k[0] in counts and k[0] in SOLVENT | IONS:
            new.append(f"{k[0]:<12s}{counts[k[0]]}")
        else:
            new.append(line)
    top.write_text("\n".join(new) + "\n")
    print(f"  wrote    {dst}  {len(written)} atoms, {len(out)} residues")


if __name__ == "__main__":
    main()
