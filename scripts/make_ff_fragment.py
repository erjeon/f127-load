#!/usr/bin/env python3
"""Collect the force field lines a guest needs and the polymer does not carry.

CHARMM-GUI returns a topology that names atom types by CGenFF name. The polymer
force field in data/toppar defines only the types the polymer uses, so a guest
that brings new ones needs them supplied alongside.

The filter is the part that has to be right, and it is wrong in both directions
if it is written carelessly. Keeping every line that mentions a new type pulls
in that type's partners, which the molecule never uses and nothing defines, so
grompp stops on the first one it reads. Keeping only the lines that mention a
new type is not enough either: the polymer force field defines fourteen bond
types in all, so a pair of ordinary types the molecule uses, CG311 with OG311
for instance, is absent as well.

So a line is kept when every type it names is one this molecule uses, and the
combination is one the polymer force field does not already define.

    python make_ff_fragment.py library/<guest> path/to/charmm36.ff [extra.itp ...]

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
from pathlib import Path

# how many atom types name an interaction in each section. Angles name three,
# not four: reading four found nothing, so the set of combinations the polymer
# already defines came out empty and every angle it carries was copied into the
# fragment as well. grompp then warned on each one.
N_NAMED = {"bondtypes": 2, "pairtypes": 2, "angletypes": 3, "dihedraltypes": 4,
           "impropertypes": 4, "cmaptypes": 5}
SECTIONS = tuple(N_NAMED)


def section(text, name):
    """Every block with this name, joined.

    A CHARMM force field writes dihedraltypes twice, once for the proper
    dihedrals and once for the impropers. Reading only the first block lost the
    two torsions doxorubicin needs around its glycosidic oxygen, and grompp then
    stopped on them with no default type.
    """
    out = []
    for piece in text.split(f"[ {name} ]")[1:]:
        out.append(piece.split("\n[")[0])
    return "\n".join(out)


def molecule_types(itp):
    out = []
    for line in itp.read_text().split("[ atoms ]")[1].split("[")[0].splitlines():
        p = line.split(";")[0].split()
        if len(p) >= 2 and p[0].isdigit():
            out.append(p[1])
    return set(out)


def defined_types(itp):
    t = itp.read_text()
    if "[ atomtypes ]" not in t:
        return set()
    return {l.split(";")[0].split()[0]
            for l in section(t, "atomtypes").splitlines() if l.split(";")[0].split()}


def key(names):
    """A combination and its reverse are the same interaction, nothing else is.

    Sorting the names instead treats every permutation as one, so a dihedral the
    polymer happened to define in another order was taken as already covered.
    Two torsions around the glycosidic oxygen of doxorubicin were dropped that
    way and grompp stopped on them.
    """
    names = tuple(names)
    return min(names, names[::-1])


def defined_combinations(itp, sec, n_named):
    """The type combinations the polymer force field already resolves."""
    t = itp.read_text()
    out = set()
    for l in section(t, sec).splitlines():
        p = l.split(";")[0].split()
        if not p:
            continue
        names = tuple(x for x in p[:n_named]
                      if not x.replace(".", "").replace("-", "").isdigit())
        if len(names) == n_named:
            out.add(key(names))
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    lib = Path(sys.argv[1])
    ff = Path(sys.argv[2])
    extras = [Path(x) for x in sys.argv[3:]]
    root = Path(__file__).resolve().parents[1]

    resn = (lib / "RESNAME.txt").read_text().strip()
    itp = next(p for p in lib.glob("*.itp") if not p.name.startswith("ff_"))
    used = molecule_types(itp)
    have = defined_types(root / "data" / "toppar" / "forcefield.itp")
    new = sorted(used - have)
    print(f"  {resn}: {len(used)} atom types, {len(new)} not carried by the polymer")
    print(f"  {' '.join(new)}")

    nb = (ff / "ffnonbonded.itp").read_text()
    lines = [f"; parameters for {resn} that data/toppar/forcefield.itp does not carry.",
             f"; new atom types: {' '.join(new)}",
             "; a line is kept only when every type it names is one this molecule",
             "; uses, so nothing is pulled in that has no definition here.",
             "", "[ atomtypes ]"]
    found = set()
    for l in section(nb, "atomtypes").splitlines():
        p = l.split(";")[0].split()
        if p and p[0] in new:
            lines.append(l.rstrip())
            found.add(p[0])
    if set(new) - found:
        sys.exit(f"  [failed] not in {ff}/ffnonbonded.itp: {sorted(set(new)-found)}")

    bodies = [(ff / "ffbonded.itp").read_text()] + [e.read_text() for e in extras]
    poly = root / "data" / "toppar" / "forcefield.itp"
    for sec in SECTIONS:
        keep, seen = [], set()
        n_named = N_NAMED[sec]
        already = defined_combinations(poly, sec, n_named)
        for body in bodies:
            for l in section(body, sec).splitlines():
                p = l.split(";")[0].split()
                if not p:
                    continue
                names = [x for x in p[:n_named] if not x.replace(".", "").replace("-", "").isdigit()]
                if not names:
                    continue
                if not all(x in used for x in names):
                    continue          # names a type this molecule never uses
                if key(names) in already:
                    continue          # the polymer force field resolves it
                if l.rstrip() in seen:
                    continue
                seen.add(l.rstrip())
                keep.append(l.rstrip())
        if keep:
            lines += ["", f"[ {sec} ]"] + keep
            print(f"  {sec:14s} {len(keep)}")

    out = lib / f"ff_{resn}.itp"
    out.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out}  {out.stat().st_size:,} B")


if __name__ == "__main__":
    main()
