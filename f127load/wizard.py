#!/usr/bin/env python3
"""Ask what system you want, then write it down so the build is reproducible.

Nothing here runs GROMACS. It collects choices, checks each one against what the
shipped structures and the geometry allow, and writes system.json. `f127 build`
reads that file, so the same system can be rebuilt later or sent to someone else
without repeating the questions.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import json
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# One micelle in a box, with a little water around it. The paper's own
# solution-route system was built at 17.0 nm, and the micelle needs 17.0 to
# stay 1.5 nm clear of its own image, so 18 leaves a margin without paying
# for water nobody reads. At 25 nm the same run costs two and a half times
# as much and answers the same question.
DEFAULT_BOX = 18.0
sys.path.insert(0, str(ROOT))
from f127load import capacity as cap
from f127load import config
from f127load import plan_size as plan
from f127load import ui

SHELL_COM_LIMIT = 1.5   # nm, set by the hollow in data/f127_shell_template.gro
MW_F127 = 12586.0
NA = 6.02214076e23
RHO = 1.0                  # g/cm3, the solution. Good to a per cent at these wt%.

# Box edge from the polymer mass fraction. Checked against both systems in the
# accompanying paper: 34 chains at 4.7 wt% gives 24.7 nm against a real 24.8, and
# at 14.5 wt% gives 17.0 nm against a real 17.0.
def box_from_wt(n_chains, wt_pct):
    m_poly = n_chains * MW_F127 / NA                 # g
    v_nm3 = (m_poly * 100.0 / wt_pct) / RHO * 1e21   # cm3 -> nm3
    return v_nm3 ** (1 / 3)


def wt_from_box(n_chains, box_nm):
    m_poly = n_chains * MW_F127 / NA
    return 100.0 * m_poly / (box_nm ** 3 * 1e-21 * RHO)


def host_extent(gro):
    """Farthest atom, the 99th percentile, and how many lie beyond it."""
    import numpy as np
    L = Path(gro).read_text().splitlines()
    n = int(L[1])
    p = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                  for l in L[2:2 + n]])
    r = np.linalg.norm(p - p.mean(0), axis=1)
    r99 = float(np.percentile(r, 99))
    return float(r.max()), r99, int((r > r99).sum())


# ------------------------------------------------------------------ prompts --
def ask(text, default=None, cast=str, check=None):
    while True:
        raw = ui.prompt(text, default)
        if not raw and default is not None:
            raw = str(default)
        if not raw:
            continue
        try:
            v = cast(raw)
        except ValueError:
            ui.note("That is not a number", "bad"); continue
        if check:
            msg = check(v)
            if msg:
                ui.note(msg, "bad"); continue
        return v


def choose(text, options, default=1):
    ui.choices(options)
    i = ask(text, default, int, lambda v: None if 1 <= v <= len(options) else "out of range")
    return options[i - 1][0]


def library():
    out = []
    for d in sorted((ROOT / "library").iterdir()):
        if not d.is_dir():
            continue
        itp = sorted(d.glob("*.itp"))
        resn = (d / "RESNAME.txt")
        crd = None
        if resn.exists():
            r = resn.read_text().strip()
            for c in (d / f"{r}.pdb", d / f"{r}.gro"):
                if c.exists():
                    crd = c
        out.append(dict(name=d.name, itp=itp[0] if itp else None,
                        crd=crd, ready=bool(itp and crd)))
    return out


def molecule_volume(entry):
    """vdW volume in nm3, from whichever structure the library has."""
    for c in (entry["crd"], *sorted((ROOT / "library" / entry["name"]).glob("*.mol2")),
              *sorted((ROOT / "library" / entry["name"]).glob("*.pdb"))):
        if c is None:
            continue
        try:
            el, xyz = cap.read_structure(Path(c))
            return cap.vdw_volume(el, xyz), cap.molecular_mass(el)
        except Exception:
            continue
    raise SystemExit(f"  [failed] could not read a structure for {entry['name']}")


# The salts the shipped topologies cover, with the anion each one brings and the
# concentration in extracellular fluid, which is where the defaults come from.
SALTS = [("NaCl", "SOD", "CLA", 1, 0.147), ("KCl", "POT", "CLA", 1, 0.0041),
         ("CaCl2", "CAL", "CLA", 2, 0.0012), ("MgCl2", "MG", "CLA", 2, 0.0009)]
PHYSIOLOGICAL = {"NaCl": 0.147, "KCl": 0.0041, "CaCl2": 0.0012, "MgCl2": 0.0009}


def ion_counts(salts, box_nm):
    """Number of each ion at the requested concentrations in this box."""
    litres = box_nm ** 3 * 1e-24
    out = {}
    for name, cat, an, z, _ in SALTS:
        c = salts.get(name, 0.0)
        if c <= 0:
            continue
        n = max(1, round(c * litres * NA))
        out[name] = dict(cation=cat, anion=an, n_cation=n, n_anion=n * z,
                         cation_charge=z)
    return out


def ask_salts(box_nm):
    """One salt, several, or the extracellular mixture.

    A single salt was all the tool allowed, so the ion series the paper actually
    ran could not be set up with it. Concentrations are asked for one salt at a
    time, because a mixture is a set of independent concentrations and not a
    single number.
    """
    how = choose("ions", [("one", "a single salt"),
                          ("mix", "a mixture, one concentration each"),
                          ("phys", "extracellular fluid, the mixture used in the paper"),
                          ("none", "none, neutralise only")], 1)
    if how == "none":
        return {}
    if how == "phys":
        salts = dict(PHYSIOLOGICAL)
    elif how == "one":
        which = choose("salt", [(s[0], s[0]) for s in SALTS], 1)
        salts = {which: ask(f"{which} concentration (M)", 0.154, float,
                            lambda v: None if 0 <= v <= 1 else "between 0 and 1")}
    else:
        salts = {}
        for name, _, _, _, default in SALTS:
            c = ask(f"{name} (M), 0 to leave it out", default, float,
                    lambda v: None if 0 <= v <= 1 else "between 0 and 1")
            if c > 0:
                salts[name] = c
    counts = ion_counts(salts, box_nm)
    if not counts:
        ui.note("no salt requested, the system will only be neutralised", "info")
        return {}
    rows = [[k, f"{salts[k]:.4f} M", f"{v['n_cation']} {v['cation']}",
             f"{v['n_anion']} {v['anion']}"] for k, v in counts.items()]
    ui.table(rows, ["salt", "concentration", "cation", "anion"])
    total = sum(v["n_cation"] + v["n_anion"] for v in counts.values())
    ui.value("ions in total", total)
    for name, v in counts.items():
        if v["n_cation"] < 5:
            ui.note(f"{name} comes to only {v['n_cation']} ion(s) in this box. A trace "
                    f"concentration needs a large box before it means anything", "warn")
    ionic = sum(salts[k] * (1 if k in ("NaCl", "KCl") else 3) for k in salts)
    ui.value("ionic strength", f"{ionic:.3f}", "M")
    return salts


def penalties(name):
    """Highest CGenFF penalty for a molecule, and the highest charge penalty.

    CHARMM-GUI writes these into the stream files it returns. They say how far
    the program had to reach by analogy: below 10 the assignment is taken from
    something closely related, 10 to 50 wants checking against something, and
    above 50 the parameter was a guess and should be fitted before the number it
    produces is trusted. The charge penalty matters more than the dihedral one
    here, because where a solute sits between a dry core and a wet corona is set
    by its partial charges.
    """
    import re
    d = ROOT / "library" / name
    worst = charge = 0.0
    for f in list(d.glob("*.rtf")) + list(d.glob("*.prm")) + list(d.glob("*.str")):
        text = f.read_text(errors="ignore")
        for m in re.finditer(r"(charge )?penalty= *([0-9.]+)", text):
            v = float(m.group(2))
            worst = max(worst, v)
            if m.group(1):
                charge = max(charge, v)
    return worst, charge


def report_penalties(name):
    worst, charge = penalties(name)
    if worst == 0:
        ui.note("no CGenFF penalties recorded. The molecule maps onto existing "
                "atom types", "ok")
        return
    tone = "bad" if max(worst, charge) > 50 else "warn" if max(worst, charge) > 10 else "good"
    ui.value("CGenFF penalty, highest", f"{worst:.1f}", tone=tone)
    ui.value("of which charge", f"{charge:.1f}", tone=tone)
    if max(worst, charge) > 50:
        ui.note("above 50 the parameter was assigned by a distant analogy. Fit it "
                "against quantum chemistry before trusting where this molecule "
                "ends up in the micelle", "bad")
    elif max(worst, charge) > 10:
        ui.note("between 10 and 50. Usable, but check the result against something "
                "measured before reporting it", "warn")


# --------------------------------------------------------------------- main --
def main():
    ui.banner("f127-load", "Design an F127 micelle system. Enter takes the value in brackets.")

    ui.step(1, 5, "The box")
    # Only the 34-chain structure exists, so the aggregation number is not a
    # choice and is no longer asked for. With the chain count fixed, the polymer
    # weight fraction is not an independent quantity either: it is the box edge
    # said another way, and offering both invited a dilute answer that tripled
    # the cost of the run for nothing. The box is what is asked, and the weight
    # fraction is reported back.
    n_chains = plan.REF_CHAINS
    g = plan.plan(n_chains)
    host = ROOT / "data" / "f127_micelle_34.gro"
    rmax = r99 = 0.0
    n_out = 0
    if host.exists():
        # The hard limit was twice the farthest atom, and that atom is the tip
        # of one stray poly(ethylene oxide) tail. Twenty-four atoms out of 69,802
        # were setting a 24 nm floor for a micelle whose bulk ends at 8.6 nm. A
        # tail that reaches past half the box folds back through the boundary,
        # which is what solvated tails do anyway, so the 99th percentile is used
        # and the few beyond it are reported rather than forbidden.
        rmax, r99, n_out = host_extent(host)
        min_box = 2 * r99 + 1.0
        built = True
    else:
        min_box, built = 2 * g["r_micelle"] + 1.0, False
    ui.value("chains", f"{n_chains}")
    ui.note("the 34-chain micelle is the only structure that ships, so the "
            "aggregation number is not asked for", "info")
    ui.value("core radius", f"{g['r_core']:.2f}", "nm")
    ui.value("micelle radius", f"{g['r_micelle']:.2f}", "nm")

    clear_box = 2 * (g["r_micelle"] + plan.MIN_CLEARANCE_NM)
    ui.note(f"one micelle sits in the box. {clear_box:.1f} nm is the smallest that "
            f"keeps it {plan.MIN_CLEARANCE_NM} nm clear of its own image, and every "
            f"nanometre past that is water", "info")
    box = ask("box edge (nm)", DEFAULT_BOX, float,
              lambda v: None if 5 <= v <= 60 else "between 5 and 60")

    # the micelle must not see its own image through the cut-off
    if box < clear_box:
        ui.note(f"that leaves {box/2 - g['r_micelle']:.2f} nm between the micelle and the "
                f"box face. {plan.MIN_CLEARANCE_NM} nm is the least that keeps it away "
                f"from its own image", "warn")
        box = clear_box
        ui.note(f"raised the box to {box:.1f} nm", "info")
    wt = wt_from_box(n_chains, box)
    ui.value("box", f"{box:.1f}", "nm")
    ui.value("clearance", f"{box/2 - g['r_micelle']:.2f}", "nm",
             "good" if box/2 - g["r_micelle"] >= plan.MIN_CLEARANCE_NM else "bad")
    ui.value("polymer", f"{wt:.1f}", "wt%")
    ui.note("the weight fraction follows from the box and is not set separately. "
            "To reach a higher one, compress a host with densify.py", "info")
    atoms = plan.plan(n_chains, box - 2 * g["r_micelle"])["atoms"]
    ui.value("atoms", f"about {atoms:,}", tone="warn" if atoms > 900_000 else None)
    if atoms > 900_000:
        ui.note(f"most of that is water. A smaller box costs less for the same "
                f"micelle, and densify.py reaches a higher concentration without "
                f"one", "warn")
    if built and box < min_box:
        ui.note(f"a {box:.1f} nm box puts the outer {n_out} atoms of the shipped "
                f"structure, which are {rmax:.1f} nm out on a few tails, past half "
                f"the box. They fold through the boundary and relax during "
                f"equilibration", "warn")
        if choose("what now", [("force", f"use {box:.1f} nm (recommended)"),
                               ("safe", f"build at {min_box:.1f} nm instead")], 1) == "safe":
            box = min_box
            wt = wt_from_box(n_chains, box)
            ui.value("box", f"{box:.1f}", "nm"); ui.value("polymer", f"{wt:.1f}", "wt%")
    if box > plan.BOX_WARN_NM:
        ui.note(f"the box is past {plan.BOX_WARN_NM:.0f} nm. This is a large system", "warn")

    ui.step(2, 5, "Ions")
    salts = ask_salts(box)

    ui.step(3, 5, "The solute")
    lib = library()
    ready = [e for e in lib if e["ready"]]
    if not ready:
        raise SystemExit("  [failed] no molecule has a topology yet. See library/README.md.")
    opts = [(e["name"], e["name"]) for e in ready]
    missing = [e["name"] for e in lib if not e["ready"]]
    guest = choose("molecule", opts, 1)
    if missing:
        ui.note(f"not ready yet: {', '.join(missing)}", "info")
    entry = next(e for e in ready if e["name"] == guest)
    v_mol, m_mol = molecule_volume(entry)
    ui.value("vdW volume", f"{v_mol:.3f}", "nm3"); ui.value("molar mass", f"{m_mol:.1f}", "g/mol")
    report_penalties(guest)

    ui.step(4, 5, "How much of it")
    mode = choose("give it as", [("count", "a number of molecules"), ("mM", "a concentration in mM"),
                                ("wt", "wt% of the polymer")], 1)
    v_core = 4 / 3 * 3.141592653589793 * g["r_core"] ** 3
    v_box_L = box ** 3 * 1e-24
    if mode == "count":
        n_guest = ask("molecules", 8, int, lambda v: None if v >= 1 else "at least 1")
    elif mode == "mM":
        c = ask("concentration (mM)", 10.0, float, lambda v: None if v > 0 else "must be positive")
        n_guest = max(1, round(c * 1e-3 * v_box_L * NA))
    else:
        w = ask("wt% of the polymer", 5.0, float, lambda v: None if v > 0 else "must be positive")
        n_guest = max(1, round(w / 100 * n_chains * MW_F127 / m_mol))
    mM = n_guest / NA / v_box_L * 1e3
    wt_guest = 100.0 * n_guest * m_mol / (n_chains * MW_F127)
    fill = 100.0 * n_guest * v_mol / v_core
    ui.value("molecules", n_guest); ui.value("concentration", f"{mM:.1f}", "mM")
    ui.value("wt% of polymer", f"{wt_guest:.1f}", "wt%")
    ui.value("core filled", f"{fill:.1f}", "%",
             "bad" if fill > 20 else "warn" if fill > 8 else "good")
    if fill > 20:
        ui.note(f"this fills {fill:.0f}% of the core. Random close packing is 64%, so it is "
                f"possible, but it is far above what is loaded experimentally", "warn")
    elif fill > 8:
        ui.note("the 100-pyrene system in the paper filled 7.4%", "info")

    ui.step(5, 5, "Where it starts")
    route = choose("start it", [("solution", "outside, dispersed in water, entering unaided"),
                            ("shell", "inside, with the shell closed around it"),
                            ("both", "both, two simulations built side by side")], 1)
    # The hollow takes far fewer molecules than the water around it, so one count
    # cannot serve both routes. Cutting the single count down to what the shell
    # holds cut the solution system with it, and the pair in the paper is 100 in
    # solution against 8 inside. The shell count is asked for on its own.
    n_shell = n_guest
    if route in ("shell", "both"):
        # The shell route places the solute inside the hollow of
        # data/f127_shell_template.gro, whose cavity is 4.6 nm across, so the
        # centre of a molecule stays within 1.5 nm.
        # The builder was asked for 54 and again for 36 and placed 32 both
        # times, which is the ceiling for pyrene in this hollow. With the
        # packing factor now set at random close packing rather than above it,
        # 0.86 of the radius gives 30, just under that ceiling, which leaves
        # room for a different random seed without being needlessly shy.
        USABLE = 0.86
        limit = SHELL_COM_LIMIT * USABLE
        n_max = 0
        while cap.cavity_radius(n_max + 1, v_mol) <= limit:
            n_max += 1
        if n_max < 1:
            ui.note(f"{guest} is too large for the hollow, so it cannot start "
                    f"inside. Using the solution route", "bad")
            route, n_shell = "solution", 0
        else:
            ui.value("the hollow holds", f"{n_max}", f"x {guest}")
            ui.note(f"the hollow is {SHELL_COM_LIMIT} nm across and random "
                    f"placement uses about {USABLE:.0%} of it", "info")
            if route == "both" and n_guest > n_max:
                ui.note(f"the solution system keeps its {n_guest}; only the "
                        f"shell system is limited", "info")
            if route == "both":
                ui.note("two simulations, not one split in half. Each is built "
                        "with its own count, and they match unless the hollow "
                        "forces otherwise", "info")
            n_shell = ask("molecules in the hollow", min(n_guest, n_max), int,
                          lambda v: None if 1 <= v <= n_max else
                          f"between 1 and {n_max}. The hollow holds no more")
            if route == "shell":
                n_guest = n_shell
            else:
                ui.value("in the water", f"{n_guest}", f"x {guest}")
                ui.value("in the hollow", f"{n_shell}", f"x {guest}")
            r_cav = cap.cavity_radius(n_shell, v_mol)
            ui.value("cavity needed", f"{r_cav:.2f}", "nm")
            ui.value("limit", f"{limit:.2f}", "nm",
                     "good" if r_cav <= limit else "bad")

    rows = [
        ("", "", "", "--- the micelle ---"),
        ("chains", n_chains, "", "the only structure that ships"),
        ("polymer", f"{wt:.2f}", "wt%",
         "follows from the box, not read. densify.py for a higher one"),
        ("box", f"{box:.2f}", "nm",
         f"{box/2 - g['r_micelle']:.2f} nm clear of the face, {plan.MIN_CLEARANCE_NM} is the least"),
        ("", "", "", ""),
        ("", "", "", "--- ions ---"),
    ]
    for name, conc in (salts or {}).items():
        c = ion_counts(salts, box).get(name, {})
        rows.append((name, f"{conc:.4f}", "M",
                     f"{c.get('n_cation','?')} {c.get('cation','')} "
                     f"and {c.get('n_anion','?')} {c.get('anion','')}"))
    if not salts:
        rows.append(("salt", "none", "", "neutralising ions only"))
    rows += [
        ("", "", "", ""),
        ("", "", "", "--- the solute ---"),
        ("guest", guest, "", f"{v_mol:.3f} nm3, {m_mol:.1f} g/mol"),
        ("count", n_guest, "",
         f"{mM:.1f} mM, {wt_guest:.1f} wt% of the polymer"),
        ("", "", "", f"fills {fill:.1f}% of the core. Experiments load 1 to 9.5 wt%"),
        ("", "", "", ""),
        ("", "", "", "--- where it starts ---"),
        ("route", route, "", "solution, shell or both"),
    ]
    if route in ("shell", "both"):
        rows.append(("inside", n_shell, "",
                     "molecules in the hollow. count is the other simulation, "
                     "not the rest of these"))
    out = Path(ask("write it to", "system.json", str))
    data = config.rows_to_json(rows, salts=tuple(ion_counts({}, 1).keys()) or
                               ("NaCl", "KCl", "CaCl2", "MgCl2"))
    config.write_json(out, data)
    cfg = {k: v for k, v in config.read(out).items() if k != "notes"}

    print()
    ui.rule()
    ui.table([[k, v] for k, v in cfg.items()])
    ui.rule()
    ui.note(f"written to {out}", "ok")
    # The config is written where the caller is standing, which is not always
    # the clone, and "./f127" only means anything inside it.
    here = Path.cwd().resolve()
    f127 = ROOT / "f127"
    cmd = "./f127" if here == ROOT else os.path.relpath(f127, here)
    if cmd.count("..") > 2:
        cmd = str(f127)
    print(f"\n    next:  {ui.bold(cmd + ' run ' + str(out))}\n")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print("\n"); ui.note("cancelled", "info")
