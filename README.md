# f127-load

Build and simulate **Pluronic F127 micelles loaded with a hydrophobic guest**, at all-atom
resolution, with GROMACS + CHARMM36/CGenFF.

This repository accompanies the manuscript *"Outer-core localization of a hydrophobic solute in Pluronic F127 micelles from all-atom molecular dynamics"*.

**[Design a system in the browser →](https://erjeon.github.io/f127-load/build.html)**
 · [What the micelle looks like](https://erjeon.github.io/f127-load/)

Both pages run entirely in the browser. Nothing is uploaded, and they work from a local
file if you would rather not use the link.

---

## What it gives you

* an **equilibrated 34-chain F127 micelle** (200 ns, CHARMM36, TIP3P) ready to load
* a **hollow-shell template** for enclosing a guest as the micelle closes
* **two loading routes**: shell closure, and spontaneous uptake from bulk water
* a **densification step** that reaches any polymer concentration without chain overlaps
* an **analysis pipeline**: radial density, encapsulation, RDF, interaction energy,
  hydration number, SASA, Rg, guest–guest association
* a **guest library** of molecules commonly loaded into F127

## What it will not do

**One aggregation number.** The 34-chain micelle is the structure that ships, and it is the
only one. Another number means a chain from CHARMM-GUI and a fresh equilibration, which is
where the time goes. The same is true of another Pluronic: the block sequence sets the
aggregation number, the core radius and the corona thickness together, so there is no
structure left to load a solute into. `plan_size.py` will tell you what building one would
cost before you start.

**One micelle per box.** The box edge is the only size you choose, and it decides how much
water surrounds the micelle. 17.0 nm is the least that keeps the micelle 1.5 nm clear of
its own periodic image; the default of 18 nm leaves a little more. Anything smaller needs
`scripts/densify.py` to compress a host first, because shrinking a box around a finished
structure folds the corona onto the core and minimisation reports an infinite force.

**Six of the twelve molecules in the library are structures only.** They are input for
parameter generation, not molecules that can be built yet. `library/README.md` says which
is which, and the tool does not offer the ones that are not ready.

It also does **not** run CGenFF itself.
[CHARMM-GUI Ligand Reader & Modeler](https://charmm-gui.org/?doc=input/ligandrm) does that
part: it looks the molecule up in the CHARMM force field library and falls back to CGenFF for
whatever is not already there, then returns the topology with a penalty on every charge and
every bonded term. Drop what it sends into `library/<name>/`. Those penalties are what
`check_params` reads. See `library/README.md`.

### What has been checked

Every molecule that has parameters was built, equilibrated and run, on one GPU:

| molecule | build | equilibrate | production | box after NpT |
|---|---|---|---|---|
| pyrene | ok | ok | ok | 17.89 nm |
| aspirin | ok | ok | ok | 17.90 nm |
| curcumin | ok | ok | ok | 17.90 nm |
| paclitaxel | ok | ok | ok | 17.90 nm |
| doxorubicin | ok | ok | ok | 17.90 nm |
| doxorubicin (+1) | ok | ok | ok | 17.90 nm |

All six were asked for an 18.0 nm box and all six settled within 0.11 nm of it.

The pyrene system was built and equilibrated a second time on a laptop, with no GPU
involved, and reached the same 17.89053 nm. That is the check that the pipeline does not
depend on the machine it was written on.

That second run is also what caught the analysis. Two of the eight analyses, the radial
distribution function and the interaction energy, had never produced anything: both name
index groups that `scripts/mkndx.py` did not write, so both stopped at `grompp` and left
their error in a log the run did not read. All eight run now, and the index carries
`core`, `corona`, `water` and `ions` alongside the groups the simulation needs.

### How long it takes

The loaded system is 577,146 atoms. Measured on that system, same box, same settings:

| | | 100 ns would take |
|---|---|---|
| one RTX 5090 | 106.3 ns/day | 0.9 days |
| one RTX 3090 | 43.1 ns/day | 2.3 days |
| one RTX 3070, 8 GB | 39.4 ns/day | 2.5 days |
| Mac mini M4, 10 CPU threads | 1.42 ns/day | 70 days |

The 5090 and 3070 rows were reported by two people who ran `./f127 run system.json 0.1`
on their own machines (2026-09-16); the other two are ours. The plan line printed at the
start of a run picks the rate for the card it finds and says so.

A laptop is for checking that the tool works, not for running it. `NPT_NS=0.02 ./f127 run
system.json 0.05` finishes the whole pipeline in under an hour and tells you that much.

---

## Install

GROMACS has to be there already; everything else is two Python packages.

```bash
git clone https://github.com/erjeon/f127-load
cd f127-load

python3 -m venv venv && source venv/bin/activate   # recent Linux and macOS
pip install "MDAnalysis>=2.8" numpy matplotlib     # refuse pip into the system Python

./f127 check                 # says what is present and what is missing
```

`./f127 check` is worth reading rather than skipping. It reports which GROMACS is on the
path and, when more than one is installed, which ones are not being used: a tpr written by
a newer GROMACS cannot be read by an older one, and finding that out from a trjconv error
takes a while. It also says whether MDAnalysis is new enough to read the tpr GROMACS 2025
writes.

## How it fits together

Clone the repository once. It is the tool, not a package you install into
somewhere else, so `./f127` is always called out of the clone. The only file you
provide is a config that names the box, the ions, the molecule and where it
starts, and everything downstream is built from it.

```
  a config                                    the clone
  system.json  ────────────►  ./f127 run  ────────────►  run_pyrene_solution/
  (anywhere on disk)                                     topology, coordinates,
                                                         tpr, trajectory, results
```

Three ways to get the config, all producing the same thing:

| | |
|---|---|
| `./f127 new` | asks the questions in the terminal |
| [the browser page](https://erjeon.github.io/f127-load/build.html) | the same choices with a sketch, then Download |
| an editor | copy the example below and change the numbers |

The config does **not** have to be copied into the clone. Give its path:

```bash
cd /path/to/f127-load
./f127 run ~/Downloads/system.json 100
```

`./f127` has to be called with a path to itself, which is what `./` means, so
running it from another directory needs the full path rather than the two
characters:

```bash
/path/to/f127-load/f127 run system.json 100
```

Either form leaves the run directory in the directory you called it from, so
results can sit on a different disk from the clone.

## Quick start

```bash
./f127 check                 # what is installed, what is missing
./f127 new                   # a few questions, writes system.json
./f127 run system.json 100   # build, equilibrate, run 100 ns
```

That is the whole sequence. `run` builds the system, minimises it, equilibrates
2 ns under NpT and then runs the production length asked for, leaving
`run_<molecule>_<route>/` behind. `./f127 build` does the same up to the tpr and
stops, for anyone who would rather submit the MD themselves.

The config is the same JSON the browser designer (`docs/build.html`) downloads,
so either source builds the same way:

```json
{
  "n_chains": 34,
  "box_nm": 18.0,
  "salts": {"NaCl": 0.147, "KCl": 0.0041},
  "guest": "pyrene",
  "n_guest": 20,
  "route": "solution"
}
```

Change a number and run `./f127 run` again. There is no need to go back through
the questions. The wizard adds a `notes` object saying where each number came
from; the build ignores it.

### Trying it without waiting a day

The default 2 ns of NpT is a few minutes on a GPU and most of a day on a laptop.
To see the whole pipeline finish first:

```bash
NPT_NS=0.02 ./f127 run system.json 0.05
```

Everything runs, nothing is equilibrated. It answers whether the tool works
here, not whether the answer is right.

### When the machine has no GPU

Nothing to set. The scripts read what `gmx --version` reports and pass the
offload flags only to a build that has them. The Homebrew GROMACS on macOS is
built with `GPU support: disabled`, and asking it for `-nb gpu` is a hard error
rather than a fallback, which is worth knowing because the message names the
flag and not the build. Thread count follows the machine as well, up to 16.

To force something different, set the flags yourself:

```bash
MDRUN_OPT="-ntmpi 1 -ntomp 8 -nb gpu -pme gpu" ./f127 run system.json 100
GMX=/usr/local/gromacs-2025.4/bin/gmx ./f127 run system.json 100
PYTHON=~/venv/bin/python ./f127 run system.json 100
```

### The other entry points

```bash
./f127 size                  # what building a different micelle would cost
./f127 pack  run_pyrene_solution   # tar.gz that rebuilds anywhere
```

### Adding a molecule

Put what CHARMM-GUI Ligand Reader and Modeler returns into `library/<name>/`, then build
the force field fragment that carries the atom types the polymer does not define:

```bash
python3 scripts/make_ff_fragment.py library/<name> /path/to/charmm36.ff
python3 scripts/check_library.py <name>      # before building anything
```

`check_library.py` answers in a second what otherwise costs a five-minute build and a
grompp error: that the residue name reads the same in all four places it appears, that the
topology and the coordinates agree atom by atom, and that every bond, angle and dihedral
the molecule needs can actually be resolved. A PDB residue name is three characters, which
is the sort of thing it catches.

Once a molecule passes, it appears in the list in `./f127 new` and in the browser page.
Molecules that have a structure but no parameters yet are left out of both, so nothing
offers you something it cannot build.

### Checking a build without running it

```bash
bash tests/smoke_test.sh pyrene     # one molecule, a few thousand steps
bash tests/edge_test.sh 3           # both ends of what the tool accepts, 3 ns each
```

The smoke test constructs the system, runs `grompp` and executes a few thousand steps,
which shows that the topology and the parameters are consistent. It does not show whether
the solute stays loaded.

The edge test takes longer and is worth the time before handing the tool to anyone else.
It builds the smallest system the tool will accept, one solute molecule with no salt in a
large box, and the largest, the shell route filled to its limit with the four-salt
physiological mixture in the smallest box that still clears the periodic image. Then it
equilibrates and runs both. Three defects that the smoke test did not reach were found
this way.

---

## Requirements

| | version | needed for |
|---|---|---|
| GROMACS | 2019 or later | building and running |
| Python | 3.8 or later | everything |
| NumPy | any recent | building and running |
| MDAnalysis | 2.0, ideally 2.8 | the analysis step only |
| matplotlib | any recent | the figures the analysis draws |

Three version notes, all of which `./f127 check` will tell you about:

**GROMACS before 2021.** The mdp files ask for the C-rescale barostat, which arrived in
2021 and is what the paper used. On an older one the run scripts substitute
Parrinello-Rahman when they set the step count, and `./f127 check` says which is in use.
Below 2019 the check stops rather than guess.

**More than one GROMACS.** The pipeline uses whichever `gmx` is on the path and never
swaps it underneath you. If a build wrote its tpr with GROMACS 2025 and an older `gmx` is
first on the path, the analysis stops with an error that does not mention versions. Source
the GMXRC of the build you mean, or set `GMX=/path/to/gmx`.

**MDAnalysis before 2.8.** GROMACS 2025 writes tpx 137, which MDAnalysis learned to read
in 2.8. On an older one the analysis reads the `.gro` beside the tpr instead and guesses
masses from atom names, which moves a centre of mass by far less than any bin width used
here. It is a note, not a failure.

Nothing else is needed. The browser pages have no dependencies: no server, no library,
no upload. They ask Google Fonts for a typeface and fall back to the system font when
there is no network, so they work offline on a cluster node.

## Citing

If you use this repository, please cite the accompanying paper and the tools it
is built on. CHARMM-GUI asks for the first three.

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
* M.J. Abraham, T. Murtola, R. Schulz, S. Páll, J.C. Smith, B. Hess, and E. Lindahl (2015).
  GROMACS: High performance molecular simulations through multi-level parallelism from
  laptops to supercomputers. *SoftwareX* 1-2:19-25. DOI: 10.1016/j.softx.2015.06.001
* N. Michaud-Agrawal, E.J. Denning, T.B. Woolf, and O. Beckstein (2011). MDAnalysis: A
  toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.*
  32:2319-2327. DOI: 10.1002/jcc.21787

## Licence

MIT.
