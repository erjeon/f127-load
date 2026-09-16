#!/usr/bin/env python3
"""Remove water that gmx solvate placed inside the solute.

solvate keeps a solvent molecule when it clears the solute by the sum of two
van der Waals radii, and it takes those radii from vdwradii.dat by atom name.
A CHARMM-GUI topology names atoms in its own way, so for anything past the
simplest molecule the radii are guessed and come out small: paclitaxel ended up
with water 0.137 nm from a heavy atom, half a bond length, and minimisation tore
the molecule apart at the first step.

Rather than raise -scale, which would thin the water everywhere including around
the polymer and change the density, the water that actually sits too close is
taken out afterwards. Whole molecules only, and the topology count follows.

    python trim_close_water.py solv.gro topol.top [cutoff_nm]

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
from pathlib import Path

import numpy as np

SOLVENT = {"TIP3", "SOL", "HOH", "WAT"}
DEFAULT_CUTOFF = 0.22          # nm, heavy atom to heavy atom


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    gro, top = Path(sys.argv[1]), Path(sys.argv[2])
    cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_CUTOFF

    lines = gro.read_text().splitlines()
    n = int(lines[1])
    head, atoms, boxline = lines[0], lines[2:2 + n], lines[2 + n]
    box = np.array([float(x) for x in boxline.split()[:3]])

    resn = np.array([a[5:10].strip() for a in atoms])
    name = np.array([a[10:15].strip() for a in atoms])
    xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])]
                    for a in atoms])
    is_wat = np.isin(resn, list(SOLVENT))
    heavy = ~np.char.startswith(np.char.upper(name.astype(str)), "H")

    solute_idx = np.where(~is_wat & heavy)[0]
    wat_idx = np.where(is_wat & heavy)[0]          # the oxygens
    if not len(solute_idx) or not len(wat_idx):
        print("  nothing to do")
        return

    # grid the solute so the search is over neighbours rather than everything
    cell = max(cutoff, 0.25)
    keys = {}
    for i in solute_idx:
        k = tuple((xyz[i] // cell).astype(int))
        keys.setdefault(k, []).append(i)

    # Which molecule each atom belongs to, by position in the file. The residue
    # number in a .gro wraps at 99,999 and a box holds 170,000 water molecules,
    # so keying on that number merged molecules that share one: 759 were found
    # and only 547 distinct keys came out, and the topology count then did not
    # match the coordinates.
    mol_of = np.empty(len(atoms), dtype=np.int64)
    m, prev = -1, None
    for i, a in enumerate(atoms):
        if a[:5] != prev:
            m += 1
            prev = a[:5]
        mol_of[i] = m

    drop_mol = set()
    for w in wat_idx:
        c = (xyz[w] // cell).astype(int)
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    near += keys.get((c[0] + dx, c[1] + dy, c[2] + dz), [])
        if not near:
            continue
        d = xyz[near] - xyz[w]
        d -= box * np.round(d / box)
        if np.linalg.norm(d, axis=1).min() < cutoff:
            drop_mol.add(int(mol_of[w]))

    if not drop_mol:
        print(f"  no water within {cutoff} nm of the solute")
        return

    keep = [a for i, a in enumerate(atoms)
            if not (a[5:10].strip() in SOLVENT and int(mol_of[i]) in drop_mol)]
    out, rid, aid, prev = [], 0, 0, None
    for a in keep:
        if a[:5] != prev:
            rid += 1
            prev = a[:5]
        aid += 1
        out.append(f"{rid % 100000:5d}{a[5:15]}{aid % 100000:5d}{a[20:]}")
    gro.write_text(f"{head}\n{len(out)}\n" + "\n".join(out) + f"\n{boxline}\n")

    removed = len(drop_mol)
    new = []
    for line in top.read_text().splitlines():
        p = line.split()
        if len(p) == 2 and p[0] in SOLVENT and p[1].isdigit():
            new.append(f"{p[0]:<12s}{int(p[1]) - removed}")
        else:
            new.append(line)
    top.write_text("\n".join(new) + "\n")
    print(f"  removed {removed} water molecule(s) closer than {cutoff} nm "
          f"to the solute")


if __name__ == "__main__":
    main()
