#!/usr/bin/env python3
"""Write the index file the run needs, straight from the .gro.

This used to go through MDAnalysis, which meant a node could not build or run a
system unless it also had the analysis stack installed. Grouping atoms by residue
name needs nothing but the fixed-format columns, so it is done here directly.
Only scripts/5_analyze.sh still needs MDAnalysis.

Atom numbers in a .gro wrap at 99,999 and a loaded micelle passes that easily, so
the numbers written here are positions in the file, counted from one, and never
the numbers printed in the file.

    python3 scripts/mkndx.py ions.gro index.ndx

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys

# These three partition the box: every atom lands in exactly one, and the check
# in main() relies on that.
GROUPS = {
    "micelle": {"PROXS", "PROXR", "PROX", "ETHOX", "ETHO"},
    "W_ION":   {"TIP3", "SOL", "SOD", "CLA", "POT", "CAL", "MG"},
}

# The two blocks separately, written in addition and overlapping the micelle, so
# they are kept out of the coverage arithmetic. scripts/5_analyze.sh names core
# as an energy group for the interaction-energy rerun and it was never written,
# so grompp stopped with "Group core referenced in the .mdp file was not found
# in the list of index groups" on every system that was ever analysed. The rerun
# behind the paper used a hand-made index that had the group.
# TIP3 was missing for the same reason: an index file supplied with -n replaces
# the moleculetype names, so "TIP3" in energygrps resolved to nothing even though
# the topology has a moleculetype by that name. Water and ions are separate here
# because an interaction energy with the solvent should not count the salt.
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
        # columns 5 to 10 hold the residue name in every gro ever written
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
