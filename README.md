# f127-load

Build and simulate a Pluronic F127 micelle loaded with a hydrophobic solute, at all-atom
resolution, with GROMACS and CHARMM36/CGenFF.

Accompanies *Outer-core localization of a hydrophobic solute in Pluronic F127 micelles from
all-atom molecular dynamics*.

[Design a system in the browser](https://erjeon.github.io/f127-load/build.html) ·
[Overview](https://erjeon.github.io/f127-load/)

## What it does

The tool ships a 34-chain PEO100–PPO65–PEO100 micelle equilibrated for 200 ns. It places the
solute either in bulk water, to watch it enter, or at the centre of a hollow PPO shell that
closes around it, then solvates, adds ions, minimises, equilibrates at 310.15 K and 1 bar and
runs production. The analysis step returns radial density, uptake counts, RDF, interaction
energies, hydration, SASA, radius of gyration and solute–solute association, with figures.

## Requirements

| | version | for |
|---|---|---|
| GROMACS | 2019 or later, GPU build | building and running |
| Python | 3.8 or later | everything |
| NumPy | any recent | building and running |
| MDAnalysis | 2.8 | the analysis step |
| matplotlib | any recent | the figures |
| CHARMM-GUI | account | only to add a molecule that is not in the library |

Throughput on the default 18 nm box (about 580,000 atoms):

| | ns/day |
|---|---|
| RTX 5090 | 106 |
| RTX 3090 | 43 |
| RTX 3070 (8 GB) | 39 |
| 10 CPU threads, no GPU | 1.4 |

## Install

```bash
git clone https://github.com/erjeon/f127-load
cd f127-load
python3 -m venv venv && source venv/bin/activate
pip install "MDAnalysis>=2.8" numpy matplotlib
./f127 check          # which GROMACS, which GPU, what is missing
```

## Run

```bash
./f127 new                        # a few questions, writes system.json
./f127 test system.json           # a short pass of every step
./f127 run  system.json 100       # build, equilibrate 2 ns, produce 100 ns
./f127 build system.json          # stop at the tpr and submit the MD yourself
./f127 pack run_pyrene_solution   # tar.gz that reruns anywhere
```

The config is the same JSON `docs/build.html` downloads, so either source builds the
same way. It can live anywhere; give its path.

```json
{
  "n_chains": 34,
  "box_nm": 18.0,
  "salts": {"NaCl": 0.147, "KCl": 0.0041},
  "solute": "pyrene",
  "n_solute": 20,
  "route": "solution"
}
```

`route` is `solution`, `shell` or `both`. Everything lands in `run_<molecule>_<route>/`:
the topology, the tpr, the xtc, and `analysis/` with the numbers and figures. Change a number
in the config and run again.

The scripts read `gmx --version` and pass GPU flags only to a build that has them. To
override:

```bash
MDRUN_OPT="-ntmpi 1 -ntomp 8 -nb gpu -pme gpu" ./f127 run system.json 100
GMX=/usr/local/gromacs-2025.4/bin/gmx ./f127 run system.json 100
NPT_NS=0.5 ./f127 run system.json 100     # equilibration length, default 2 ns
```

## Molecules

Ready to run: pyrene, doxorubicin (neutral and protonated), paclitaxel, curcumin, aspirin.
Structure only, parameters still to be made: docetaxel, ibuprofen, indomethacin, nile red,
quercetin, resveratrol. `python3 scripts/check_library.py` lists the state of every entry.

To add one, run [CHARMM-GUI Ligand Reader and Modeler](https://charmm-gui.org/?doc=input/ligandrm),
put the result into `library/<name>/`, then

```bash
python3 scripts/make_ff_fragment.py library/<name> /path/to/charmm36.ff
python3 scripts/check_library.py <name>
python3 -m f127load.check_params library/<name>   # CGenFF penalties
bash tests/smoke_test.sh <name>
```

`library/README.md` has the file list and the steps.

## Tests

```bash
bash tests/smoke_test.sh pyrene     # build, grompp, a few thousand steps
bash tests/edge_test.sh 3           # smallest and largest system the tool accepts, 3 ns each
```

## References

* S. Jo, T. Kim, V.G. Iyer, and W. Im (2008). CHARMM-GUI: A Web-based Graphical User
  Interface for CHARMM. *J. Comput. Chem.* 29:1859-1865. DOI: 10.1002/jcc.20945
* J. Lee, X. Cheng, J.M. Swails, M.S. Yeom, P.K. Eastman, J.A. Lemkul, S. Wei, J. Buckner,
  J.C. Jeong, Y. Qi, S. Jo, V.S. Pande, D.A. Case, C.L. Brooks III, A.D. MacKerell Jr,
  J.B. Klauda, and W. Im (2016). CHARMM-GUI Input Generator for NAMD, GROMACS, AMBER,
  OpenMM, and CHARMM/OpenMM Simulations using the CHARMM36 Additive Force Field.
  *J. Chem. Theory Comput.* 12:405-413. DOI: 10.1021/acs.jctc.5b00935
* S. Kim, J. Lee, S. Jo, C.L. Brooks III, H.S. Lee, and W. Im (2017). CHARMM-GUI Ligand
  Reader and Modeler for CHARMM Force Field Generation of Small Molecules.
  *J. Comput. Chem.* 38:1879-1886. DOI: 10.1002/jcc.24829
* K. Vanommeslaeghe, E. Hatcher, C. Acharya, S. Kundu, S. Zhong, J. Shim, E. Darian,
  O. Guvench, P. Lopes, I. Vorobyov, and A.D. MacKerell Jr (2010). CHARMM General Force
  Field: A force field for drug-like molecules compatible with the CHARMM all-atom additive
  biological force fields. *J. Comput. Chem.* 31:671-690. DOI: 10.1002/jcc.21367
* K. Vanommeslaeghe and A.D. MacKerell Jr (2012). Automation of the CHARMM General Force
  Field (CGenFF) I: bond perception and atom typing. *J. Chem. Inf. Model.* 52:3144-3154.
  DOI: 10.1021/ci300363c
* K. Vanommeslaeghe, E.P. Raman, and A.D. MacKerell Jr (2012). Automation of the CHARMM
  General Force Field (CGenFF) II: assignment of bonded parameters and partial atomic
  charges. *J. Chem. Inf. Model.* 52:3155-3168. DOI: 10.1021/ci3003649
* M.J. Abraham, T. Murtola, R. Schulz, S. Páll, J.C. Smith, B. Hess, and E. Lindahl (2015).
  GROMACS: High performance molecular simulations through multi-level parallelism from
  laptops to supercomputers. *SoftwareX* 1-2:19-25. DOI: 10.1016/j.softx.2015.06.001
* N. Michaud-Agrawal, E.J. Denning, T.B. Woolf, and O. Beckstein (2011). MDAnalysis: A
  toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.*
  32:2319-2327. DOI: 10.1002/jcc.21787

## Licence

MIT.
