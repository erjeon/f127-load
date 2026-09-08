#!/usr/bin/env python3
"""Report how far a CGenFF parameter set had to guess.

CGenFF assigns a penalty to every parameter it could not take directly from
its training set. A low penalty means an exact or near analogue was available.
A high one means the value was extrapolated and may be wrong.

The usual reading is that below 10 the parameter can be used as it stands,
between 10 and 50 it deserves a look, and above 50 it should be checked against
a quantum calculation before any quantitative claim rests on it.

Where the penalty sits matters as much as how large it is. A guessed torsion
changes how freely a group rotates. A guessed charge changes how the molecule
sees water, and for a question about where a solute sits in a micelle that is
the more serious of the two.

  python3 check_params.py path/to/charmm-gui-output/

Copyright (c) Pukyong National University / NCHM Lab
Eunryul Jeon  <qlsguswjs@pukyong.ac.kr>
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

SECTIONS = ("BONDS", "ANGLES", "DIHEDRALS", "IMPROPERS", "NONBONDED")
GOOD, WATCH = 10.0, 50.0


def collect(root: Path):
    """Return {kind: [(penalty, text)]} over every prm and rtf under root."""
    found = defaultdict(list)
    for path in list(root.rglob("*.prm")) + list(root.rglob("*.rtf")):
        # the CHARMM-GUI toppar directory holds the stock force field, which
        # carries no penalties of its own and would only add noise
        if "toppar" in path.parts:
            continue
        kind = None
        for line in path.read_text(errors="ignore").splitlines():
            head = line.strip().split(" ")[0].upper()
            if head in SECTIONS:
                kind = head
                continue
            # the RESI line of the rtf carries the two summary figures for the
            # whole molecule, the worst bonded penalty and the worst charge one
            summary = re.search(
                r"param penalty=\s*([\d.]+).*?charge penalty=\s*([\d.]+)", line)
            if summary:
                found["_summary"] = (float(summary.group(1)),
                                     float(summary.group(2)))
                continue
            m = re.search(r"penalty\s*=\s*([\d.]+)", line)
            if not m:
                continue
            value = float(m.group(1))
            if line.lstrip().upper().startswith("ATOM"):
                found["CHARGES"].append((value, line.split("!")[0].strip()))
            elif kind:
                found[kind].append((value, line.split("!")[0].strip()))
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="f127load-check-params",
        description="Report CGenFF penalties for a parameter set.",
        epilog="Pukyong National University / NCHM Lab.  "
               "Eunryul Jeon <qlsguswjs@pukyong.ac.kr>")
    ap.add_argument("directory", type=Path,
                    help="CHARMM-GUI output directory for one molecule")
    ap.add_argument("--list", type=float, default=None, metavar="P",
                    help="also list every term whose penalty exceeds P")
    args = ap.parse_args(argv)

    if not args.directory.exists():
        sys.exit(f"not found: {args.directory}")
    found = collect(args.directory)
    summary = found.pop("_summary", None)
    if not found and not summary:
        sys.exit("no penalties found. Is this a CHARMM-GUI output directory?")

    print(f"{args.directory.name}\n")
    print(f"  {'term':<12}{'count':>7}{'median':>9}{'max':>8}"
          f"{'>10':>6}{'>50':>6}")
    worst_charge = summary[1] if summary else 0.0
    total_bad = 0
    for kind in ("CHARGES", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPERS"):
        vals = sorted(v for v, _ in found.get(kind, []))
        if not vals:
            continue
        over = sum(1 for v in vals if v > WATCH)
        total_bad += over
        if kind == "CHARGES":
            worst_charge = max(worst_charge, vals[-1])
        print(f"  {kind.lower():<12}{len(vals):>7}{vals[len(vals)//2]:>9.1f}"
              f"{vals[-1]:>8.1f}{sum(1 for v in vals if v > GOOD):>6}{over:>6}")

    if summary:
        print(f"\n  CHARMM-GUI summary for the whole molecule")
        print(f"    worst bonded penalty {summary[0]:.1f}")
        print(f"    worst charge penalty {summary[1]:.1f}")
        if summary[1] > WATCH:
            total_bad += 1

    print()
    if total_bad == 0:
        print("  Every term was assigned from a close analogue. No action needed.")
    else:
        print(f"  {total_bad} term(s) above {WATCH:.0f}. These were extrapolated.")
        if worst_charge > WATCH:
            print(f"  The charges reach {worst_charge:.0f}, which is the one to worry")
            print("  about. Partial charges set how the molecule sees water, so a")
            print("  guessed charge moves the balance between polar and non-polar")
            print("  and can move where the solute ends up sitting.")
        else:
            print("  The charges are within range, so the electrostatics are on")
            print("  firmer ground than the bonded terms. Guessed torsions change")
            print("  how freely a group rotates, which matters less for a question")
            print("  about partitioning than for one about conformation.")
        print()
        print("  Reasonable responses, cheapest first.")
        print("    report the penalties and keep the molecule to qualitative use")
        print("    check the partition coefficient against experiment")
        print("    refit the offending terms against a quantum calculation")

    if args.list is not None:
        print(f"\n  terms above {args.list:g}")
        rows = [(v, k, t) for k, lst in found.items() for v, t in lst
                if v > args.list]
        for v, k, t in sorted(rows, key=lambda r: -r[0]):
            print(f"    {v:6.1f}  {k.lower():<10} {t[:56]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
