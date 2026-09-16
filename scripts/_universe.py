"""Open a trajectory, falling back to a .gro when MDAnalysis cannot read the tpr.

GROMACS 2025 writes tpx 137, which MDAnalysis reads only from 2.8. The analyses
here need only residue names, atom names and masses, so a .gro from the same run
stands in and masses are guessed from atom names.
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
