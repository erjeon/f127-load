#!/usr/bin/env python3
"""Write the index file the run needs, straight from the .gro.

Grouping atoms by residue name needs nothing but the fixed-format columns, so
no MDAnalysis is needed here.

Atom numbers in a .gro wrap at 99,999 and a loaded micelle passes that easily, so
the numbers written here are positions in the file, counted from one, and never
the numbers printed in the file.

    python3 scripts/mkndx.py ions.gro index.ndx
"""
import sys

# These three partition the box: every atom lands in exactly one, and the check
# in main() relies on that.
GROUPS = {
    "micelle": {"PROXS", "PROXR", "PROX", "ETHOX", "ETHO"},
    "W_ION":   {"TIP3", "SOL", "SOD", "CLA", "POT", "CAL", "MG"},
}

# Written in addition to the groups above and overlapping them, so they are
# kept out of the coverage check. scripts/5_analyze.sh uses core and water as
# energy groups; ions are separate so a solvent interaction energy has no salt.
BLOCKS = {
    "core":   {"PROXS", "PROXR", "PROX"},
    "corona": {"ETHOX", "ETHO"},
    "water":  {"TIP3", "SOL"},
    "ions":   {"SOD", "CLA", "POT", "CAL", "MG"},
}


def read_resnames(path):
    with open(path) as fh:
        fh.readline()                       # title
        n = int(fh.readline())
        # columns 5 to 10 hold the residue name
        return [fh.readline()[5:10].strip() for _ in range(n)]


def write_group(fh, name, members):
    fh.write(f"[ {name} ]\n")
    for i, a in enumerate(members):
        fh.write(f"{a:6d}" + ("\n" if (i + 1) % 15 == 0 else " "))
    if len(members) % 15:
        fh.write("\n")


def main():
    gro, out = sys.argv[1], sys.argv[2]
    res = read_resnames(gro)

    # whatever is left over after the micelle and the solvent is the solute, so a
    # new molecule needs no edit here
    named = GROUPS["micelle"] | GROUPS["W_ION"]
    solute = sorted({r for r in res if r not in named})
    if len(solute) > 1:
        sys.exit(f"[failed] more than one candidate solute residue: {solute}")

    groups = {"micelle": GROUPS["micelle"], "W_ION": GROUPS["W_ION"]}
    if solute:
        groups[solute[0]] = {solute[0]}

    picked = {k: [i + 1 for i, r in enumerate(res) if r in v] for k, v in groups.items()}
    total = sum(len(v) for v in picked.values())
    if total != len(res):
        missing = sorted({r for r in res} - named - set(solute))
        sys.exit(f"[failed] groups cover {total} of {len(res)} atoms, unassigned: {missing}")

    blocks = {k: [i + 1 for i, r in enumerate(res) if r in v]
              for k, v in BLOCKS.items()}

    with open(out, "w") as fh:
        write_group(fh, "System", range(1, len(res) + 1))
        for k, v in picked.items():
            write_group(fh, k, v)
        for k, v in blocks.items():
            write_group(fh, k, v)
    print("  ndx ok:", {k: len(v) for k, v in {**picked, **blocks}.items()})


if __name__ == "__main__":
    main()
