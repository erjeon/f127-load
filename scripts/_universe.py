"""Open a trajectory whose topology MDAnalysis may not be able to read.

GROMACS 2025 writes tpx 137. MDAnalysis only learned to read that in 2.8, and
2.4 ships with several distributions still, so the analysis half of this tool
died on its own output with

    NotImplementedError: Your tpx version is 137, which this parser does not
    support, yet

A run directory always holds a .gro beside the .tpr, and the analyses here need
only residue names, atom names and masses. So the tpr is tried first, for its
exact masses, and a .gro stands in when it cannot be read. Masses are then
guessed from atom names, which moves a centre of mass by far less than the bin
width used anywhere in this tool.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
import warnings
from pathlib import Path

import MDAnalysis as mda


def open_universe(topology, trajectory):
    """Universe(topology, trajectory), falling back to a .gro in the same run."""
    top = Path(topology)
    tried = []
    candidates = [top]
    if top.suffix == ".tpr":
        stem = top.with_suffix("")
        candidates += [p for p in (top.parent / "npt.gro", top.parent / "prod.gro",
                                   top.parent / "em.gro", stem.with_suffix(".gro"))
                       if p.exists()]
    for cand in candidates:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(cand), str(trajectory))
            if cand != top:
                print(f"  [note] {top.name} could not be read by MDAnalysis "
                      f"{mda.__version__}; using {cand.name} for the topology. "
                      f"Masses are guessed from atom names.")
            return u
        except Exception as exc:
            tried.append(f"{cand.name}: {type(exc).__name__}")
    sys.exit("  [failed] no readable topology. Tried " + "; ".join(tried) +
             f".\n  MDAnalysis {mda.__version__} reads tpx 137 only from 2.8 "
             f"onward. Upgrade it, or keep npt.gro in the run directory.")
