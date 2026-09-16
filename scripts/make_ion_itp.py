#!/usr/bin/env python3
"""Turn the CHARMM ion parameters into the GROMACS files the tool ships.

Only sodium and chloride came with the original build, so the ion series that
the paper actually ran, potassium and the divalent cations, could not be set up
with the tool. The parameters are in the CHARMM stream file that CHARMM-GUI
supplies with every job, and they only need converting.

CHARMM writes epsilon in kcal/mol as a negative number and Rmin/2 in angstrom.
GROMACS wants sigma in nm and a positive epsilon in kJ/mol, with

    sigma = 2 * (Rmin/2) / 2^(1/6) / 10        Rmin/2 in A -> sigma in nm
    epsilon = |eps_kcal| * 4.184

The NBFIX pairs matter as much as the types. CHARMM corrects several cation
oxygen pairs away from the combination rule, and dropping them changes how
tightly a cation sits on a poly(ethylene oxide) ether oxygen, which is the
quantity the ion part of the paper is about.

    python make_ion_itp.py <toppar_water_ions.str> <outdir> [forcefield.itp]

A pair is only written if both of its types are in the shipped force field.
grompp refuses a nonbond_params line that names a type the system does not
define, and for F127 none of the corrections apply anyway: they are written for
carboxylate, ester and phosphate oxygens, and the ether oxygen of poly(ethylene
oxide) is OG301, which appears in none of them. Cation and ether oxygen
therefore interact through the plain combination rule here.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import re
import sys
from pathlib import Path

KCAL = 4.184
SIG = 2.0 / 2.0 ** (1 / 6) / 10.0          # (Rmin/2 in A) -> sigma in nm

IONS = {                                    # resname: (charge, mass, element)
    "SOD": (1.0, 22.98977, "Na"), "POT": (1.0, 39.09830, "K"),
    "CLA": (-1.0, 35.45000, "Cl"), "CAL": (2.0, 40.08000, "Ca"),
    "MG": (2.0, 24.30500, "Mg"),
}
ATNUM = {"Na": 11, "K": 19, "Cl": 17, "Ca": 20, "Mg": 12}


def parse(path):
    text = Path(path).read_text(errors="ignore")
    nb, fix, sect = {}, [], None
    for line in text.splitlines():
        s = line.split("!")[0].rstrip()
        if not s.strip():
            continue
        head = s.split()[0].upper()
        if head in ("NONBONDED", "NBFIX", "BONDS", "ANGLES", "DIHEDRALS",
                    "IMPROPER", "CMAP", "HBOND", "END", "READ", "ATOMS"):
            sect = head
            continue
        f = s.split()
        if sect == "NONBONDED" and len(f) >= 4 and f[0].upper() in IONS:
            nb[f[0].upper()] = (abs(float(f[2])) * KCAL, float(f[3]) * SIG)
        elif sect == "NBFIX" and len(f) >= 4:
            a, b = f[0].upper(), f[1].upper()
            if a in IONS or b in IONS:
                fix.append((a, b, abs(float(f[2])) * KCAL, float(f[3]) * SIG / 2))
    return nb, fix


def have_types(path):
    if not path:
        return None
    keep, on = set(), False
    for line in Path(path).read_text().splitlines():
        s = line.split(";")[0].strip()
        if s.startswith("["):
            on = s.replace(" ", "").startswith("[atomtypes]")
            continue
        if on and s and len(s.split()) > 5:
            keep.add(s.split()[0])
    return keep


def main():
    src, out = sys.argv[1], Path(sys.argv[2])
    ff = sys.argv[3] if len(sys.argv) > 3 else None
    nb, fix = parse(src)
    known = have_types(ff)
    if known is not None:
        known |= set(IONS)
        dropped = [(a, b) for a, b, _, _ in fix if a not in known or b not in known]
        fix = [f for f in fix if f[0] in known and f[1] in known]
        if dropped:
            miss = sorted({t for pair in dropped for t in pair} - known)
            print(f"  dropped {len(dropped)} NBFIX pairs naming types this system does "
                  f"not have: {', '.join(miss)}")
    out.mkdir(parents=True, exist_ok=True)
    missing = [k for k in IONS if k not in nb]
    if missing:
        raise SystemExit(f"  [failed] no parameters for {', '.join(missing)}")

    # whatever the force field already defines must not be defined again, or
    # grompp stops on a duplicate atom type
    already = set()
    pairs_already = set()
    if ff:
        sect = None
        for line in Path(ff).read_text().splitlines():
            s = line.split(";")[0].strip()
            if s.startswith("["):
                sect = s.replace(" ", "")
                continue
            f = s.split()
            if sect == "[atomtypes]" and len(f) > 5:
                already.add(f[0])
            elif sect == "[nonbond_params]" and len(f) >= 4:
                pairs_already.add(frozenset((f[0], f[1])))
        nb = {k: v for k, v in nb.items() if k not in already}
        fix = [f for f in fix if frozenset((f[0], f[1])) not in pairs_already]
        if already & set(IONS):
            print(f"  already in the force field, not repeated: "
                  f"{', '.join(sorted(already & set(IONS)))}")

    types = ["[ atomtypes ]",
             "; name  at.num   mass     charge  ptype     sigma       epsilon"]
    for name, (eps, sigma) in sorted(nb.items()):
        q, m, el = IONS[name]
        types.append(f"{name:>8} {ATNUM[el]:5d} {m:10.5f} {q:10.4f}   A "
                     f"{sigma:18.12e} {eps:14.6e}")
    pairs = ["", "[ nonbond_params ]", "; i        j    func      sigma        epsilon"]
    for a, b, eps, sigma in fix:
        pairs.append(f"{a:>8} {b:>7}     1 {sigma:18.12e} {eps:14.6e}")
    (out / "ion_types.itp").write_text("\n".join(types + pairs) + "\n")

    for name, (q, m, el) in IONS.items():
        if (out / f"{name}.itp").exists() and name in already:
            continue
        (out / f"{name}.itp").write_text(
            f"; {name}, converted from CHARMM by make_ion_itp.py\n\n"
            "[ moleculetype ]\n; name  nrexcl\n"
            f"{name}       1\n\n"
            "[ atoms ]\n; nr  type  resnr  residu  atom  cgnr  charge      mass\n"
            f"    1  {name:>5}      1  {name:>6}  {name:>5}     1  {q:7.3f}  {m:10.5f}\n")
    print(f"  {len(nb)} ion types, {len(fix)} NBFIX pairs -> {out}")
    for name, (eps, sigma) in sorted(nb.items()):
        print(f"    {name:4s} sigma {sigma:.5f} nm  epsilon {eps:.4f} kJ/mol")


if __name__ == "__main__":
    main()
