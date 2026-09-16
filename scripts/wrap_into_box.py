#!/usr/bin/env python3
"""Fold whole molecules back inside the box, and say how many atoms moved.

A structure equilibrated in a large box keeps a few chain tips two or three nm
outside a smaller one. gmx editconf resizes and centres but leaves them there,
gmx solvate then puts water where they will be once the periodic image is taken,
and minimisation starts with an infinite force.

Molecules are moved as a whole, by the same shift for every atom in them, so no
bond is stretched across the boundary. Molecule boundaries come from the residue
numbering in the file, which for this pipeline is one number per chain.

    python wrap_into_box.py in.gro out.gro

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
from pathlib import Path

import numpy as np


def main():
    src, dst = sys.argv[1], sys.argv[2]
    L = Path(src).read_text().splitlines()
    n = int(L[1])
    rows = L[2:2 + n]
    box = np.array([float(x) for x in L[2 + n].split()[:3]])
    pos = np.array([[float(r[20:28]), float(r[28:36]), float(r[36:44])] for r in rows])
    resid = np.array([int(r[0:5]) for r in rows])

    moved = 0
    out = pos.copy()
    for r in np.unique(resid):
        m = resid == r
        com = pos[m].mean(0)
        shift = -box * np.floor(com / box)
        if np.any(shift):
            out[m] = pos[m] + shift
            moved += int(m.sum())

    lines = [f"{r[:20]}{out[i,0]:8.3f}{out[i,1]:8.3f}{out[i,2]:8.3f}"
             for i, r in enumerate(rows)]
    Path(dst).write_text(L[0] + "\n" + L[1] + "\n" + "\n".join(lines) + "\n"
                         + L[2 + n] + "\n")
    print(moved)


if __name__ == "__main__":
    main()
