#!/usr/bin/env python3
"""Estimate how many copies of a solute a Pluronic F127 micelle can hold.

The micelle core is not an empty container.  Every solute molecule placed
inside has to displace poly(propylene oxide), so the number that fits is set
by how much of the core volume may be given up before the assembly stops
behaving like a micelle.  Three separate numbers are reported.

  recommended   loading at 5 wt% of the polymer mass, the range most often
                reported for F127 formulations
  high          loading at 10 wt%, at the upper end of what is reported
  ceiling       the geometric limit at which the cavity needed to hold the
                solutes reaches 75% of the core radius.  Above this the
                shell cannot close around them and the run will not
                equilibrate.  This is a hard stop, not a recommendation.

Volumes are van der Waals volumes computed by Monte Carlo integration over
the union of atomic spheres using Bondi radii.

Copyright (c) Pukyong National University / NCHM Lab
Eunryul Jeon  <qlsguswjs@pukyong.ac.kr>
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

# Bondi van der Waals radii, in angstrom
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
}

ATOMIC_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904,
}

# Reference F127 micelle shipped with this tool
N_CHAINS = 34
MW_F127 = 12586.0            # g/mol, PEO100–PPO65–PEO100
R_CORE_NM = 3.99             # PPO 90th percentile radius, measured on
                             # data/f127_micelle_34.gro. The accompanying paper
                             # quotes 3.91 nm for a different system, the one
                             # loaded from solution in a 17 nm box.
PO_PER_CHAIN = 65
MW_PO = 58.08                # g/mol, one propylene oxide unit
RHO_PPO = 1.004              # g/cm3, bulk poly(propylene oxide)
NA = 6.02214076e23

# A solute cluster whose enclosing sphere exceeds this fraction of the core
# radius cannot be wrapped by the shell during relaxation.
CAVITY_LIMIT = 0.75

# A cavity has to be larger than the molecules it holds. 1.4 assumes they end up
# at 71 per cent of the volume, which is between random close packing at 64 and
# the crystalline limit at 74, so it describes molecules that were stacked
# rather than dropped. place_guest drops them at random and rejects overlaps,
# and asked for 54 pyrene in the shell hollow it placed 32. 1.0/0.64 is what
# random close packing actually allows.
PACKING_FACTOR = 1.56


def read_mol2(path: Path):
    """Return (elements, coordinates in angstrom) from a TRIPOS mol2 file."""
    lines = path.read_text().splitlines()
    try:
        start = lines.index("@<TRIPOS>ATOM") + 1
    except ValueError:
        raise ValueError(f"{path} has no @<TRIPOS>ATOM record")
    stop = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("@<TRIPOS>"):
            stop = i
            break
    elements, coords = [], []
    for line in lines[start:stop]:
        fields = line.split()
        if len(fields) < 6:
            continue
        elements.append(fields[5].split(".")[0])
        coords.append([float(fields[2]), float(fields[3]), float(fields[4])])
    if not elements:
        raise ValueError(f"{path} contains no atom records")
    return elements, np.asarray(coords)


def read_pdb(path: Path):
    """Return (elements, coordinates in angstrom) from a PDB file."""
    elements, coords = [], []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = line[76:78].strip()
        if not element:
            element = line[12:16].strip()[0]
        elements.append(element.capitalize() if len(element) > 1 else element.upper())
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not elements:
        raise ValueError(f"{path} contains no atom records")
    return elements, np.asarray(coords)


def read_structure(path: Path):
    if path.suffix.lower() == ".mol2":
        return read_mol2(path)
    if path.suffix.lower() == ".pdb":
        return read_pdb(path)
    raise ValueError(f"unsupported structure format: {path.suffix}")


def vdw_volume(elements, coords, n_samples=400_000, seed=0):
    """Van der Waals volume in nm3, by Monte Carlo over the union of spheres."""
    radii = np.array([VDW_RADII.get(e, 1.70) for e in elements])
    lo = (coords - radii[:, None]).min(axis=0)
    hi = (coords + radii[:, None]).max(axis=0)
    box = hi - lo
    rng = np.random.default_rng(seed)
    points = lo + rng.random((n_samples, 3)) * box
    inside = np.zeros(n_samples, dtype=bool)
    # chunk over atoms to keep the distance array small
    for i in range(0, len(radii), 64):
        c = coords[i:i + 64]
        r = radii[i:i + 64]
        d2 = ((points[:, None, :] - c[None, :, :]) ** 2).sum(axis=2)
        inside |= (d2 <= (r ** 2)[None, :]).any(axis=1)
    volume_a3 = box.prod() * inside.mean()
    return volume_a3 / 1000.0          # angstrom^3 to nm^3


def molecular_mass(elements):
    return sum(ATOMIC_MASS.get(e, 12.011) for e in elements)


def gyration_radius(coords):
    centre = coords.mean(axis=0)
    return float(np.sqrt(((coords - centre) ** 2).sum(axis=1).mean())) / 10.0


def bounding_radius(coords, elements):
    """Radius in nm of the smallest sphere about the centroid holding the molecule."""
    radii = np.array([VDW_RADII.get(e, 1.70) for e in elements])
    centre = coords.mean(axis=0)
    d = np.sqrt(((coords - centre) ** 2).sum(axis=1)) + radii
    return float(d.max()) / 10.0


def cavity_radius(n_copies, molecule_volume_nm3):
    """Radius of the sphere needed to hold n copies packed together."""
    total = n_copies * molecule_volume_nm3 * PACKING_FACTOR
    return (3.0 * total / (4.0 * math.pi)) ** (1.0 / 3.0)


def core_properties():
    polymer_mass_g = N_CHAINS * MW_F127 / NA
    ppo_mass_g = N_CHAINS * PO_PER_CHAIN * MW_PO / NA
    ppo_volume_nm3 = ppo_mass_g / RHO_PPO * 1e21
    geometric_volume_nm3 = 4.0 / 3.0 * math.pi * R_CORE_NM ** 3
    return polymer_mass_g, ppo_volume_nm3, geometric_volume_nm3


def copies_at_weight_fraction(fraction, solute_mw):
    """Number of copies giving `fraction` of the polymer mass."""
    polymer_mass_g, _, _ = core_properties()
    solute_mass_g = solute_mw / NA
    return int(round(polymer_mass_g * fraction / solute_mass_g))


def geometric_ceiling(molecule_volume_nm3, single_radius_nm):
    """Largest copy count whose packed cavity stays inside CAVITY_LIMIT * R_core."""
    limit = CAVITY_LIMIT * R_CORE_NM
    if single_radius_nm > limit:
        return 0
    n = 1
    while cavity_radius(n + 1, molecule_volume_nm3) <= limit:
        n += 1
        if n > 100_000:
            break
    return n


def analyse(path: Path, quiet=False):
    elements, coords = read_structure(path)
    volume = vdw_volume(elements, coords)
    mw = molecular_mass(elements)
    r_single = bounding_radius(coords, elements)
    rec = copies_at_weight_fraction(0.05, mw)
    high = copies_at_weight_fraction(0.10, mw)
    ceiling = geometric_ceiling(volume, r_single)
    polymer_mass_g, ppo_volume, geom_volume = core_properties()
    return {
        "name": path.stem,
        "n_atoms": len(elements),
        "mw": mw,
        "volume_nm3": volume,
        "radius_nm": r_single,
        "rg_nm": gyration_radius(coords),
        "recommended": min(rec, ceiling),
        "high": min(high, ceiling),
        "ceiling": ceiling,
        "rec_capped": rec > ceiling,
        "high_capped": high > ceiling,
        "fill_at_rec": min(rec, ceiling) * volume / ppo_volume,
        "wt_at_ceiling": ceiling * (mw / NA) / polymer_mass_g,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="f127load-capacity",
        description="Estimate how many copies of a solute an F127 micelle can hold.",
        epilog="Pukyong National University / NCHM Lab.  Eunryul Jeon "
               "<qlsguswjs@pukyong.ac.kr>",
    )
    parser.add_argument("structures", nargs="+", type=Path,
                        help="solute structure files (.mol2 or .pdb)")
    parser.add_argument("--csv", type=Path, default=None,
                        help="also write the table to this CSV file")
    args = parser.parse_args(argv)

    polymer_mass_g, ppo_volume, geom_volume = core_properties()
    rows = []
    for path in args.structures:
        if not path.exists():
            print(f"[skip] {path} not found", file=sys.stderr)
            continue
        try:
            rows.append(analyse(path))
        except ValueError as exc:
            print(f"[skip] {path}: {exc}", file=sys.stderr)

    rows.sort(key=lambda r: r["mw"])

    print(f"Reference micelle : {N_CHAINS} F127 chains, R_core = {R_CORE_NM} nm")
    print(f"PPO block volume  : {ppo_volume:.1f} nm3 "
          f"(geometric core sphere {geom_volume:.1f} nm3)")
    print(f"Polymer mass      : {N_CHAINS * MW_F127 / 1000:.1f} kg/mol per micelle")
    print()
    header = (f"{'solute':<14}{'MW':>8}{'vol':>8}{'r':>7}"
              f"{'rec':>6}{'high':>6}{'ceil':>7}{'fill':>7}")
    print(header)
    print(f"{'':<14}{'g/mol':>8}{'nm3':>8}{'nm':>7}{'5wt%':>6}{'10wt%':>6}{'max':>7}{'%':>7}")
    print("-" * len(header))
    for r in rows:
        flag = "*" if r["rec_capped"] else " "
        print(f"{r['name']:<14}{r['mw']:>8.1f}{r['volume_nm3']:>8.3f}"
              f"{r['radius_nm']:>7.2f}{r['recommended']:>6d}{r['high']:>6d}"
              f"{r['ceiling']:>7d}{100 * r['fill_at_rec']:>6.1f}{flag}")
    print()
    print("rec   copies at 5 wt% of polymer mass, the usual experimental range")
    print("high  copies at 10 wt%, upper end of reported loadings")
    print("ceil  geometric limit, cavity reaches 75% of the core radius")
    print("fill  fraction of the PPO block volume filled at the recommended count")
    print("*     the weight-fraction value was cut down to the geometric ceiling")

    if args.csv:
        import csv
        with args.csv.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["solute", "n_atoms", "mw_g_per_mol", "vdw_volume_nm3",
                             "radius_nm", "rg_nm", "copies_5wt", "copies_10wt",
                             "copies_ceiling", "core_fill_fraction_at_5wt",
                             "wt_fraction_at_ceiling"])
            for r in rows:
                writer.writerow([r["name"], r["n_atoms"], f"{r['mw']:.2f}",
                                 f"{r['volume_nm3']:.4f}", f"{r['radius_nm']:.3f}",
                                 f"{r['rg_nm']:.3f}", r["recommended"], r["high"],
                                 r["ceiling"], f"{r['fill_at_rec']:.4f}",
                                 f"{r['wt_at_ceiling']:.4f}"])
        print(f"\nwritten: {args.csv.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
