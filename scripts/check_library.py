#!/usr/bin/env python3
"""Check every molecule in the library before anything is built.

Each of these went wrong at least once while doxorubicin was being added, and
each cost a build of several minutes before grompp said so. They are all
answerable from the files themselves, in a second.

    python scripts/check_library.py [name ...]

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POLY = ROOT / "data" / "toppar" / "forcefield.itp"
N_NAMED = {"bondtypes": 2, "angletypes": 3, "dihedraltypes": 4}
MOL_SEC = {"bonds": ("bondtypes", 2), "angles": ("angletypes", 3),
           "dihedrals": ("dihedraltypes", 4)}


def sections(text, name):
    """Every block with this name. CHARMM writes dihedraltypes twice."""
    return "\n".join(p.split("\n[")[0] for p in text.split(f"[ {name} ]")[1:])


def rows(text, name):
    out = []
    for line in sections(text, name).splitlines():
        p = line.split(";")[0].split()
        if p:
            out.append(p)
    return out


def key(names):
    names = tuple(names)
    return min(names, names[::-1])


def combos(text, sec, n):
    out = set()
    for p in rows(text, sec):
        names = tuple(x for x in p[:n]
                      if not x.replace(".", "").replace("-", "").isdigit())
        if len(names) == n:
            out.add(key(names))
    return out


def check(lib):
    """(name, problems, resn, atoms, charge). problems is None for a molecule
    that was never prepared, which is a different thing from a broken one: half
    the library ships as a structure and nothing else, and printing six
    [failed] lines for that reads as damage rather than as work not done."""
    name = lib.name
    problems = []
    resfile = lib / "RESNAME.txt"
    itps = [q for q in lib.glob("*.itp") if not q.name.startswith("ff_")]
    if not resfile.exists() and not itps:
        return name, None, "-", 0, 0.0
    if not resfile.exists():
        return name, ["no RESNAME.txt"]
    resn = resfile.read_text().strip()
    if len(resn) > 3:
        problems.append(f"RESNAME {resn} is longer than a PDB residue name allows")

    if not itps:
        return name, ["no topology"]
    itp = itps[0]
    text = itp.read_text()

    mt = rows(text, "moleculetype")
    declared = mt[0][0] if mt else "?"
    if declared != resn:
        problems.append(f"moleculetype is {declared}, RESNAME.txt says {resn}")

    atoms = rows(text, "atoms")
    itp_res = {a[3] for a in atoms if len(a) > 3}
    if itp_res != {resn}:
        problems.append(f"topology labels residues {sorted(itp_res)}")
    names = [a[4] for a in atoms if len(a) > 4]
    used = {a[1] for a in atoms if len(a) > 1}
    charge = sum(float(a[6]) for a in atoms if len(a) > 6)

    crd = next((p for p in (lib / f"{resn}.pdb", lib / f"{resn}.gro",
                            lib / f"{name}.pdb", lib / f"{name}.gro")
                if p.exists()), None)
    if crd is None:
        problems.append(f"no {resn}.pdb")
    else:
        if crd.suffix == ".pdb":
            lines = [l for l in crd.read_text().splitlines()
                     if l.startswith(("ATOM", "HETATM"))]
            cnames = [l[12:16].strip() for l in lines]
            cres = {l[17:20].strip() for l in lines}
        else:
            L = crd.read_text().splitlines()
            n = int(L[1])
            cnames = [l[10:15].strip() for l in L[2:2 + n]]
            cres = {l[5:10].strip() for l in L[2:2 + n]}
        if len(cnames) != len(names):
            problems.append(f"{crd.name} has {len(cnames)} atoms, topology has {len(names)}")
        elif cnames != names:
            i = next(i for i, (a, b) in enumerate(zip(names, cnames)) if a != b)
            problems.append(f"atom {i+1} is {names[i]} in the topology, {cnames[i]} in {crd.name}")
        if cres != {resn}:
            problems.append(f"{crd.name} labels residues {sorted(cres)}")

    # can every interaction be resolved from the polymer plus the fragment
    frag = lib / f"ff_{resn}.itp"
    have_t = POLY.read_text()
    if frag.exists():
        have_t += "\n" + frag.read_text()
    defined = {p[0] for p in rows(have_t, "atomtypes")}
    missing_types = sorted(used - defined)
    if missing_types:
        problems.append(f"atom types with no definition: {' '.join(missing_types)}")
    else:
        ty = {a[0]: a[1] for a in atoms if len(a) > 1}
        for mol_sec, (ff_sec, n) in MOL_SEC.items():
            have = combos(have_t, ff_sec, n)
            miss = set()
            for p in rows(text, mol_sec):
                if len(p) < n + 1 or not p[0].isdigit() or len(p) > n + 2:
                    continue
                try:
                    k = key(tuple(ty[x] for x in p[:n]))
                except KeyError:
                    continue
                if k not in have:
                    miss.add(k)
            if miss:
                problems.append(f"{len(miss)} {mol_sec} with no parameters, "
                                f"first {' '.join(sorted(miss)[0])}")
    return name, problems, resn, len(names), charge


def main():
    wanted = sys.argv[1:]
    libs = sorted(p for p in (ROOT / "library").iterdir()
                  if p.is_dir() and (not wanted or p.name in wanted))
    bad = ready = plain = 0
    print(f"{'molecule':<20}{'resn':<6}{'atoms':>6}{'charge':>8}  state")
    for lib in libs:
        out = check(lib)
        name, problems = out[0], out[1]
        extra = out[2:] if len(out) > 2 else ("?", 0, 0.0)
        resn, natoms, charge = extra if len(extra) == 3 else ("?", 0, 0.0)
        # a neutral molecule sums to something like -4e-9, and "-0.00" reads
        # as a charge that is not quite zero rather than as rounding
        if abs(charge) < 5e-3:
            charge = 0.0
        num = f"{natoms:>6}{charge:>8.2f}" if problems is not None else f"{'':>14}"
        if problems is None:
            plain += 1
            print(f"{name:<20}{resn:<6}{num}  structure only, no parameters yet")
        elif problems:
            bad += 1
            print(f"{name:<20}{resn:<6}{num}  [failed]")
            for p in problems:
                print(f"{'':34}{p}")
        else:
            ready += 1
            print(f"{name:<20}{resn:<6}{num}  ok")
    print(f"\n  {ready} of {len(libs)} ready to build", end="")
    print(f", {plain} structure only" if plain else "", end="")
    print(f", {bad} broken" if bad else "")
    if plain and not wanted:
        print("  A structure-only molecule needs a CHARMM-GUI job. "
              "See library/README.md")
    # Only a half-prepared molecule is a failure. Six of these shipped as
    # structures on purpose, so they must not make the check exit non-zero.
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
