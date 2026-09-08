#!/usr/bin/env python3
"""Place solute molecules in or around a Pluronic F127 micelle.

Three loading modes are available.

  soak   The solutes are dispersed in the aqueous phase, clear of the micelle,
         and are left to enter on their own during the simulation.  Nothing is
         assumed about where they end up.  This is the unbiased route.

  core   A cavity is opened at the centre of the micelle and the solutes are
         packed into it.  The cavity is sized from the van der Waals volume of
         the solutes, so it is no larger than it has to be.  Polymer monomers
         inside the cavity are displaced outward as rigid units rather than
         deleted, which keeps every chain intact, and the shell is closed back
         onto the solutes by the staged relaxation in 3_relax.sh.

  both   The core is filled first, the placement is verified, and any remaining
         copies are then dispersed in the aqueous phase.  Use this to study
         exchange between an already loaded core and free solute in solution.

The number of copies that a micelle can reasonably hold is not obvious.  Run
capacity.py first, which reports a recommended count at 5 wt% of the polymer
mass, a high count at 10 wt%, and the geometric ceiling above which the shell
cannot close.

Copyright (c) Pukyong National University / NCHM Lab
Eunryul Jeon  <qlsguswjs@pukyong.ac.kr>
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
except ImportError:  # pragma: no cover
    sys.exit(
        "MDAnalysis is required and was not found.\n"
        "    python3 -m venv .venv && source .venv/bin/activate\n"
        '    pip install "MDAnalysis>=2.0"\n'
        "Recent Linux and macOS installs refuse pip into the system Python "
        "(PEP 668), so use a virtual environment."
    )

from capacity import (CAVITY_LIMIT, PACKING_FACTOR, VDW_RADII, cavity_radius,
                      read_structure, vdw_volume)

PPO_RESNAMES = ("PROX", "PROXR", "PROXS")
PEO_RESNAMES = ("ETHO", "ETHOX")
WATER_RESNAMES = ("TIP3", "SOL", "HOH", "WAT")

CLASH_CUTOFF = 2.2       # angstrom, minimum heavy atom separation on placement
SOAK_MARGIN = 15.0       # angstrom, closest approach of a soaked solute to the micelle
CAVITY_PAD = 3.0         # angstrom, extra room around the packed solutes


def log(message):
    print(f"[load] {message}")


def core_radius(ppo_positions, centre, percentile=90.0):
    d = np.linalg.norm(ppo_positions - centre, axis=1)
    return float(np.percentile(d, percentile))


def random_rotation(rng):
    """Uniform random rotation matrix by QR of a Gaussian matrix."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def has_clash(candidate, occupied, box, cutoff=CLASH_CUTOFF):
    if len(occupied) == 0:
        return False
    pairs = capped_distance(candidate, occupied, max_cutoff=cutoff,
                            box=box, return_distances=False)
    return len(pairs) > 0


def place_in_cavity(template, n_copies, centre, r_cavity, box, rng, occupied):
    """Pack n_copies of `template` inside a sphere of radius r_cavity."""
    placed = []
    template = template - template.mean(axis=0)
    reach = np.linalg.norm(template, axis=1).max()
    usable = max(r_cavity - reach, 0.0)
    for index in range(n_copies):
        for attempt in range(4000):
            # sample uniformly in the usable sphere, biased to the centre early on
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            radius = usable * rng.random() ** (1.0 / 3.0)
            coords = template @ random_rotation(rng).T + centre + direction * radius
            pool = occupied if len(placed) == 0 else np.vstack([occupied] + placed)
            if not has_clash(coords, pool, box):
                placed.append(coords)
                break
        else:
            raise RuntimeError(
                f"could not place copy {index + 1} of {n_copies} inside the cavity. "
                "Reduce --n, or run capacity.py to see the ceiling for this solute."
            )
    return placed


def place_in_solution(template, n_copies, centre, r_exclude, box, rng, occupied):
    """Scatter n_copies of `template` in the aqueous phase outside r_exclude."""
    placed = []
    template = template - template.mean(axis=0)
    reach = np.linalg.norm(template, axis=1).max()
    half = box[:3] / 2.0
    for index in range(n_copies):
        for attempt in range(8000):
            point = rng.random(3) * box[:3]
            if np.linalg.norm(point - centre) < r_exclude + reach:
                continue
            # keep the whole molecule inside the box so nothing straddles an edge
            if np.any(point - reach < 0) or np.any(point + reach > box[:3]):
                continue
            coords = template @ random_rotation(rng).T + point
            pool = occupied if len(placed) == 0 else np.vstack([occupied] + placed)
            if not has_clash(coords, pool, box):
                placed.append(coords)
                break
        else:
            raise RuntimeError(
                f"could not place copy {index + 1} of {n_copies} in solution. "
                "The box may be too small or too crowded."
            )
    return placed


def displace_polymer(universe, centre, r_cavity):
    """Move polymer monomers out of the cavity as rigid units.

    Atoms are not deleted.  Every residue whose centre falls inside the cavity
    is translated radially outward until its centre sits just outside, which
    leaves the internal geometry of each monomer untouched and only stretches
    the bonds between neighbouring monomers.  Energy minimisation recovers
    those.  Deleting atoms instead would sever the chains.
    """
    polymer = universe.select_atoms(
        " or ".join(f"resname {name}" for name in PPO_RESNAMES + PEO_RESNAMES))
    if len(polymer) == 0:
        raise RuntimeError("no polymer residues found. Check the residue names "
                           "in the micelle file.")
    moved = 0
    for residue in polymer.residues:
        positions = residue.atoms.positions
        centroid = positions.mean(axis=0)
        offset = centroid - centre
        distance = np.linalg.norm(offset)
        if distance >= r_cavity:
            continue
        direction = offset / distance if distance > 1e-6 else np.array([0.0, 0.0, 1.0])
        residue.atoms.positions = positions + direction * (r_cavity - distance)
        moved += 1
    return moved, len(polymer.residues)


def water_to_remove(universe, centre, r_cavity):
    """Water residues whose oxygen sits inside the cavity."""
    water = universe.select_atoms(
        " or ".join(f"resname {name}" for name in WATER_RESNAMES))
    if len(water) == 0:
        return water
    d = np.linalg.norm(water.positions - centre, axis=1)
    inside = water[d < r_cavity]
    return inside.residues.atoms


def build_argument_parser():
    parser = argparse.ArgumentParser(
        prog="f127load-place",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Place solute molecules in or around an F127 micelle.",
        epilog="""\
examples
  # unbiased: 100 copies dispersed in water, let them find their own place
  python3 load.py --micelle data/f127_micelle_34.gro \\
                  --solute library/pyrene/pyrene.pdb --mode soak --n 100

  # loaded core: 58 copies packed inside, shell closed by 3_relax.sh
  python3 load.py --micelle data/f127_micelle_34.gro \\
                  --solute library/curcumin/curcumin.pdb --mode core --n 58

  # both: 25 inside and 25 free in solution, for exchange studies
  python3 load.py --micelle data/f127_micelle_34.gro \\
                  --solute library/paclitaxel/paclitaxel.pdb \\
                  --mode both --n 50 --n-core 25

Run capacity.py first to see how many copies are reasonable for your solute.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
""")
    parser.add_argument("--micelle", type=Path, required=True,
                        help="equilibrated micelle structure (.gro or .pdb)")
    parser.add_argument("--solute", type=Path, required=True,
                        help="single solute structure (.pdb or .mol2)")
    parser.add_argument("--mode", choices=("soak", "core", "both"), required=True,
                        help="where to put the solutes")
    parser.add_argument("--n", type=int, required=True,
                        help="total number of copies")
    parser.add_argument("--n-core", type=int, default=None,
                        help="copies to put inside when --mode both "
                             "(default: half of --n)")
    parser.add_argument("--resname", default="LIG",
                        help="residue name written for the solute (default LIG)")
    parser.add_argument("--out", type=Path, default=Path("loaded.gro"),
                        help="output structure (default loaded.gro)")
    parser.add_argument("--seed", type=int, default=20260903,
                        help="random seed, so a placement can be reproduced")
    parser.add_argument("--force", action="store_true",
                        help="place more copies than the geometric ceiling allows. "
                             "The shell will probably fail to close.")
    return parser


def main(argv=None):
    args = build_argument_parser().parse_args(argv)
    rng = np.random.default_rng(args.seed)

    universe = mda.Universe(str(args.micelle))
    box = universe.dimensions
    if box is None or np.any(box[:3] <= 0):
        sys.exit("the micelle file has no box. Add one with gmx editconf.")

    ppo = universe.select_atoms(
        " or ".join(f"resname {name}" for name in PPO_RESNAMES))
    if len(ppo) == 0:
        sys.exit("no PPO residues found. Expected one of: "
                 + ", ".join(PPO_RESNAMES))
    centre = ppo.center_of_mass()
    r_core = core_radius(ppo.positions, centre)
    log(f"micelle: {len(universe.atoms)} atoms, box {box[0] / 10:.2f} nm")
    log(f"core   : centre {centre / 10} nm, R_core {r_core / 10:.2f} nm")

    elements, template = read_structure(args.solute)
    volume = vdw_volume(elements, template)
    log(f"solute : {args.solute.stem}, {len(elements)} atoms, "
        f"vdW volume {volume:.3f} nm3")

    n_core = 0
    n_soak = 0
    if args.mode == "core":
        n_core = args.n
    elif args.mode == "soak":
        n_soak = args.n
    else:
        n_core = args.n_core if args.n_core is not None else args.n // 2
        n_soak = args.n - n_core
        if n_core < 0 or n_soak < 0:
            sys.exit("--n-core must be between 0 and --n")
    log(f"mode   : {args.mode}  ({n_core} inside, {n_soak} in solution)")

    solute_atoms = []
    removed = universe.atoms[[]]

    if n_core > 0:
        r_cavity = cavity_radius(n_core, volume) * 10.0 + CAVITY_PAD
        fraction = r_cavity / r_core
        log(f"cavity : radius {r_cavity / 10:.2f} nm, "
            f"{100 * fraction:.0f}% of R_core")
        if fraction > CAVITY_LIMIT:
            message = (f"the cavity needed for {n_core} copies is "
                       f"{100 * fraction:.0f}% of the core radius, above the "
                       f"{100 * CAVITY_LIMIT:.0f}% limit. The shell will not "
                       f"close. Run capacity.py to see the ceiling.")
            if not args.force:
                sys.exit("[FAIL] " + message)
            log("WARNING (forced): " + message)

        removed = water_to_remove(universe, centre, r_cavity)
        if len(removed):
            log(f"water  : removing {len(removed.residues)} residues "
                f"from the cavity")
        moved, total = displace_polymer(universe, centre, r_cavity)
        log(f"polymer: displaced {moved} of {total} monomers out of the cavity "
            f"(none deleted)")

        keep = universe.atoms.difference(removed)
        solute_atoms += place_in_cavity(template, n_core, centre, r_cavity,
                                        box, rng, keep.positions)
        log(f"placed : {n_core} copies inside the cavity")

    if n_soak > 0:
        keep = universe.atoms.difference(removed)
        occupied = keep.positions
        if solute_atoms:
            occupied = np.vstack([occupied] + solute_atoms)
        exclude = r_core + SOAK_MARGIN
        solute_atoms += place_in_solution(template, n_soak, centre, exclude,
                                          box, rng, occupied)
        log(f"placed : {n_soak} copies in solution, at least "
            f"{SOAK_MARGIN / 10:.1f} nm clear of the core")

    write_output(universe, removed, solute_atoms, elements, args, box)
    return 0


def write_output(universe, removed, solute_atoms, elements, args, box):
    keep = universe.atoms.difference(removed)
    n_solute_atoms = len(elements) * len(solute_atoms)
    total = len(keep) + n_solute_atoms

    names = []
    for index, element in enumerate(elements, start=1):
        names.append(f"{element}{index}"[:5])

    lines = [f"F127 micelle loaded with {len(solute_atoms)} x "
             f"{args.solute.stem} ({args.mode} mode)", f"{total:5d}"]

    # existing atoms keep their identity
    for atom in keep:
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (
            atom.resid % 100000, atom.resname[:5], atom.name[:5],
            atom.index % 100000 + 1,
            atom.position[0] / 10, atom.position[1] / 10, atom.position[2] / 10))

    resid = (keep.residues.resids.max() if len(keep) else 0)
    serial = len(keep)
    for coords in solute_atoms:
        resid += 1
        for name, position in zip(names, coords):
            serial += 1
            lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (
                resid % 100000, args.resname[:5], name, serial % 100000,
                position[0] / 10, position[1] / 10, position[2] / 10))

    lines.append("%10.5f%10.5f%10.5f" % (box[0] / 10, box[1] / 10, box[2] / 10))
    args.out.write_text("\n".join(lines) + "\n")

    log(f"wrote  : {args.out.resolve()}  ({total} atoms)")

    note = args.out.with_suffix(".topology.txt")
    water_removed = len(removed.residues) if len(removed) else 0
    note.write_text(
        "Update the [ molecules ] section of your topology to match this "
        "structure.\n\n"
        f"  add     {args.resname:<8} {len(solute_atoms)}\n"
        f"  subtract water residues removed from the cavity: {water_removed}\n\n"
        "The solute itp must be included before [ system ]. Parameters for a "
        "new molecule can be generated with CHARMM-GUI Ligand Reader and "
        "Modeler, which returns CGenFF parameters compatible with CHARMM36.\n"
    )
    log(f"wrote  : {note.resolve()}")


if __name__ == "__main__":
    sys.exit(main())
