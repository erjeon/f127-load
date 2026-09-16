#!/usr/bin/env python3
"""Print the settings the build script needs, one line, space separated.

Kept out of 1_build.sh because a here-document inside a process substitution is
awkward to get right and impossible to test on its own. This can be run by hand:

    python scripts/read_config.py system.json

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from f127load import config

NA = 6.02214076e23
SALTS = {"NaCl": ("SOD", "CLA", 1), "KCl": ("POT", "CLA", 1),
         "CaCl2": ("CAL", "CLA", 2), "MgCl2": ("MG", "CLA", 2)}


def main():
    # A mistyped path used to come back as a nine-line FileNotFoundError
    # traceback out of pathlib, which reads as a broken tool rather than as a
    # name that does not exist.
    if len(sys.argv) < 2:
        sys.exit("usage: read_config.py system.json")
    if not Path(sys.argv[1]).is_file():
        sys.exit(f"[failed] no such config: {sys.argv[1]}")
    c = config.read(sys.argv[1])
    chains = c.get("chains", c.get("n_chains"))
    box = float(c.get("box", c.get("box_nm")))
    guest = c.get("guest")
    n_guest = c.get("count", c.get("n_guest"))
    route = c.get("route")

    counts = c.get("ion_counts")
    if not counts:
        counts = {}
        litres = box ** 3 * 1e-24
        # A molarity reaches here two ways. system.cfg writes it at the top level,
        # one line per salt; the browser page and anything written by hand put it
        # inside a "salts" object. Only the first was read, so a config carrying
        # {"salts": {"NaCl": 0.154}} built a system with no salt in it and said
        # nothing, and genion's "No ions to add" went to a log.
        molar = {}
        nested = c.get("salts")
        if isinstance(nested, dict):
            molar.update(nested)
        for k in SALTS:
            if isinstance(c.get(k), (int, float)):
                molar[k] = c[k]
        unknown = sorted(set(molar) - set(SALTS))
        if unknown:
            print(f"[warn] {sys.argv[1]}: no ion parameters for {', '.join(unknown)}. "
                  f"Known salts are {', '.join(SALTS)}", file=sys.stderr)
        for k, (cat, an, z) in SALTS.items():
            v = molar.get(k)
            if isinstance(v, (int, float)) and v > 0:
                n = max(1, round(v * litres * NA))
                counts[k] = dict(cation=cat, anion=an, n_cation=n,
                                 n_anion=n * z, cation_charge=z)
        if not counts and isinstance(c.get("salt_M"), (int, float)) \
                and c["salt_M"] > 0:                       # an older JSON config
            n = max(1, round(c["salt_M"] * litres * NA))
            counts["NaCl"] = dict(cation="SOD", anion="CLA", n_cation=n,
                                  n_anion=n, cation_charge=1)
        if not counts and molar and not unknown:
            print(f"[warn] {sys.argv[1]} asks for salt but every value is zero or "
                  "not a number, so the system will be built with none",
                  file=sys.stderr)
    spec = ",".join(f"{k}:{v['cation']}:{v['anion']}:{v['n_cation']}:"
                    f"{v['n_anion']}:{v.get('cation_charge', 1)}"
                    for k, v in counts.items()) or "-"
    first = next(iter(counts), "none")
    conc = 0.0
    if first != "none":
        nested = c.get("salts") if isinstance(c.get("salts"), dict) else {}
        conc = c.get(first, nested.get(first, c.get("salt_M", 0.0)))
    # How many start inside the hollow, which the shell route holds far fewer of
    # than the water around it. A settings file written before this existed says
    # nothing, and then both systems take the same count, which is what used to
    # happen.
    n_inside = c.get("inside", n_guest)
    print(chains, box, first, conc, guest, n_guest, route, spec, n_inside)


if __name__ == "__main__":
    main()
