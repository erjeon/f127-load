# f127-load

Build and simulate **Pluronic F127 micelles loaded with a hydrophobic guest**, at all-atom
resolution, with GROMACS + CHARMM36/CGenFF.

This repository accompanies the manuscript *"Entry of a hydrophobic solute into the core
interior of Pluronic F127 micelles is slow: an all-atom molecular dynamics study"*.

---

## What it gives you

* an **equilibrated 34-chain F127 micelle** (200 ns, CHARMM36, TIP3P) ready to load
* a **hollow-shell template** for enclosing a guest as the micelle closes
* **two loading routes** — shell closure and spontaneous uptake from bulk water
* a **densification step** that reaches any polymer concentration without chain overlaps
* an **analysis pipeline** — radial density, encapsulation, RDF, interaction energy,
  hydration number, SASA, Rg, guest–guest association
* a **guest library** of molecules commonly loaded into F127

## What it does not do

It does **not** generate CGenFF parameters. Use
[CHARMM-GUI Ligand Reader & Modeler](https://charmm-gui.org/?doc=input/ligandrm) and drop the
downloaded topology into `library/<name>/`. See `library/README.md`.

---

## Quick start

```bash
# 0. put the guest topology in place (once, from CHARMM-GUI)
#    library/ibuprofen/IBU.itp   +   RESNAME.txt containing "IBU"

# 1. build a loaded micelle at 4.7 wt%
bash scripts/2_load.sh --guest ibuprofen --n 8 --method shell --box 24.8 --out run_ibu

# 2. equilibrate
bash scripts/3_equilibrate.sh run_ibu

# 3. (optional) densify to a higher polymer concentration
python3 scripts/densify.py run_ibu/npt1.gro run_ibu/c15.gro run_ibu/topol.top 0.15

# 4. production
bash scripts/4_production.sh run_ibu 100      # ns

# 5. analyse + plot
bash scripts/5_analyze.sh run_ibu --begin 50000
```

A build can be checked without running any production MD:

```bash
bash tests/smoke_test.sh ibuprofen
```

This constructs the system, runs `grompp`, and executes a few thousand MD steps to confirm
the topology and parameters are consistent. It does **not** test whether the guest stays
loaded.

---

## Requirements

GROMACS 2022 or later, Python 3.8+ with NumPy and MDAnalysis, Open Babel (for 3D structure
generation from SMILES only).

## Citing

If you use this repository, please cite the accompanying paper and:

* S. Kim, *et al.*, CHARMM-GUI ligand reader and modeler for CHARMM force field generation of
  small molecules, *J. Comput. Chem.* **38**, 1879 (2017). DOI: 10.1002/jcc.24829
* K. Vanommeslaeghe, *et al.*, CHARMM general force field, *J. Comput. Chem.* **31**, 671
  (2010). DOI: 10.1002/jcc.21367
* M. J. Abraham, *et al.*, GROMACS, *SoftwareX* **1–2**, 19 (2015).
  DOI: 10.1016/j.softx.2015.06.001

## Licence

MIT.
