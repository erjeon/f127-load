#!/usr/bin/env python3
"""How large a micelle of N chains ends up being.

The aggregation number of F127 is not a single number. It rises with temperature
and differs between preparations, so anyone reusing this tool will want a
micelle other than the 34 chains that ship with it. This module gives the
geometry and the system size before the box is built.

Every relation here is either geometry or is calibrated against a system that
was actually run, and the calibration is stated where it is used. No value is
taken from the literature: the experimental range of the aggregation number is
for the user to supply.
"""
import argparse

# The micelle distributed with this tool, data/f127_micelle_34.gro
REF_CHAINS = 34
REF_R_CORE = 3.99          # nm, PPO 90th percentile, measured on that file
REF_R_MICELLE = 7.0        # nm, where the corona runs out
ATOMS_PER_CHAIN = 2053

# Calibrated on a 400 ns solution-loaded run: 34 chains in a 17.0 nm box came to
# 488,823 atoms, of which 419,021 water and ions at 33.4 water/nm3 fill 4,182 of
# the 4,913 nm3, so the polymer excludes 731 nm3.
ATOMS_PER_NM3_WATER = 100.1
EXCLUDED_NM3_PER_REF = 731.0

# Past this the system is mostly water. It is a note, not a limit.
BOX_WARN_NM = 22.0

# Clearance between the micelle surface and the box face. Below this the
# micelle sees its own periodic image through the 1.2 nm cut-off.
MIN_CLEARANCE_NM = 1.5
WATER_PAD_NM = 2 * MIN_CLEARANCE_NM


def plan(n_chains, pad_nm=WATER_PAD_NM):
    """Geometry and system size for a micelle of n_chains."""
    scale = (n_chains / REF_CHAINS) ** (1 / 3)      # core volume goes as N
    r_core = REF_R_CORE * scale
    # The corona is one PEO block long however many chains there are, so its
    # thickness is carried across rather than scaled.
    r_micelle = r_core + (REF_R_MICELLE - REF_R_CORE)
    box = 2 * r_micelle + pad_nm

    polymer = ATOMS_PER_CHAIN * n_chains
    excluded = EXCLUDED_NM3_PER_REF * n_chains / REF_CHAINS
    water = ATOMS_PER_NM3_WATER * max(box ** 3 - excluded, 0.0)
    atoms = int(round(polymer + water))

    return dict(n_chains=n_chains, r_core=r_core, r_micelle=r_micelle, box=box,
                atoms=atoms, polymer=polymer, large=box > BOX_WARN_NM)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("n_chains", nargs="*", type=int, default=[20, 34, 50, 75, 100],
                    help="aggregation numbers to compare (default: 20 34 50 75 100)")
    ap.add_argument("--pad", type=float, default=WATER_PAD_NM,
                    help=f"water across the gap to the periodic image, nm "
                         f"(default {WATER_PAD_NM})")
    a = ap.parse_args()

    print(f"{'chains':>7}{'R core':>9}{'R micelle':>11}{'box':>8}{'atoms':>11}")
    print(f"{'':>7}{'nm':>9}{'nm':>11}{'nm':>8}")
    print("-" * 46)
    for n in a.n_chains:
        p = plan(n, a.pad)
        mark = "  <- over %.0f nm" % BOX_WARN_NM if p["large"] else ""
        star = " *" if n == REF_CHAINS else "  "
        print(f"{n:>7}{p['r_core']:>9.2f}{p['r_micelle']:>11.2f}{p['box']:>8.1f}"
              f"{p['atoms']:>11,}{star}{mark}")
    print("\n* the micelle that ships with this tool, data/f127_micelle_34.gro")
    print("Any other aggregation number has to be built by the shell-closure route,")
    print("scripts/2_load.sh, with the chain count changed. The template scales, the")
    print("equilibration does not: give a larger micelle proportionally longer to relax.")
    print(f"\nA box past {BOX_WARN_NM:.0f} nm is a large system.")
    print("\nThe aggregation number itself is an input, not a prediction. It rises with")
    print("temperature and differs between preparations, so take it from a measurement")
    print("on the batch being modelled.")


if __name__ == "__main__":
    main()
