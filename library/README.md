# Guest library

Molecules that are commonly loaded into Pluronic F127 micelles, prepared for use with this
workflow. Properties are taken from PubChem (CID given); XLogP is the PubChem computed value.

3D structures (`*.mol2`, `*.pdb`) were generated from the canonical SMILES with Open Babel
(`--gen3d --minimize --ff MMFF94`) and are provided **only as input for parameter
generation**. They are not equilibrated.

## Status

| molecule | CID | formula | MW | XLogP | atoms | charge | topology | tier |
|---|---|---|---|---|---|---|---|---|
| pyrene            | 31423   | C16H10      | 202.3 | 4.9 | 26  |  0 | `PYR.itp` | A |
| aspirin           | 2244    | C9H8O4      | 180.2 | 1.2 | 21  |  0 | `ASP.itp` | A |
| curcumin          | 969516  | C21H20O6    | 368.4 | 3.2 | 47  |  0 | `CUR.itp` | A |
| paclitaxel        | 36314   | C47H51NO14  | 853.9 | 2.5 | 113 |  0 | `TAX.itp` | B |
| doxorubicin       | 31703   | C27H29NO11  | 543.5 | 1.3 | 68  |  0 | `DOX.itp` | C |
| doxorubicin (+1)  | 31703   | C27H30NO11+ | 544.5 | 1.3 | 69  | +1 | `DXP.itp` | C |
| ibuprofen         | 3672    | C13H18O2    | 206.3 | 3.5 | 33  |  0 | structure only | A |
| indomethacin      | 3715    | C19H16ClNO4 | 357.8 | 4.3 | 41  |  0 | structure only | A |
| nile_red          | 65182   | C20H18N2O2  | 318.4 | 3.8 | 42  |  0 | structure only | A |
| resveratrol       | 445154  | C14H12O3    | 228.2 | 3.1 | 29  |  0 | structure only | A |
| quercetin         | 5280343 | C15H10O7    | 302.2 | 1.5 | 32  |  0 | structure only | B |
| docetaxel         | 148124  | C43H53NO14  | 807.9 | 1.6 | 111 |  0 | structure only | B |

The six with a topology are the six the tool offers. The rest carry a structure
and nothing else, so `./f127 new` and the browser page leave them out rather
than let a build fail five minutes in. `python3 scripts/check_library.py` prints
this state for every molecule in a second.

**Tier A** is small, neutral and drug-like, and CGenFF penalties are expected to
be low. **Tier B** is large or highly substituted, so the penalty report is worth
reading; some dihedrals may need attention. **Tier C** is charged.

Doxorubicin ships twice on purpose. `DOX` is the neutral molecule as CHARMM-GUI
returns it. `DXP` carries the protonated daunosamine amine, which is the species
at physiological pH, and its net charge is +1.0000 rather than something close
to it. The loading mechanism in the accompanying paper is dispersive, so a
charged guest is expected to behave differently; that is the reason to have both
rather than a reason to avoid one.

Note that the carboxylic acids (ibuprofen, aspirin, indomethacin) are supplied as the
**neutral** forms. At pH 7.4 they are largely deprotonated; the neutral form is the species
that partitions into the core, but this should be stated explicitly in any work using them.

## Adding a topology

1. Open [CHARMM-GUI Ligand Reader & Modeler](https://charmm-gui.org/?doc=input/ligandrm).
2. Upload `library/<name>/<name>.mol2` (or paste the SMILES from `SMILES.txt`).
3. Select **GROMACS** output. Download the job archive.
4. From the archive, copy into `library/<name>/`:
   * `toppar/<RESI>.itp`, the guest topology
   * `ligandrm.pdb` as `<RESI>.pdb`, the coordinates from the same job
   * the CGenFF penalty report, saved as `penalty.txt`
5. Create `RESNAME.txt` containing the residue name used in the itp, `IBU` say.
   Three characters. A PDB residue field is three columns wide and a fourth
   spills into the chain identifier, which is how `DOXP` became `DXP`.
6. Build the force field fragment, then check the whole thing:

   ```bash
   python3 scripts/make_ff_fragment.py library/<name> /path/to/charmm36.ff
   python3 scripts/check_library.py <name>
   ```

7. Record the CHARMM-GUI job id in `CHARMMGUI_JOBID.txt` for reproducibility.
8. Verify with `bash tests/smoke_test.sh <name>`.

**Check the penalty report.** CGenFF reports a penalty for every assigned parameter. Values
below 10 are generally acceptable; above 50 the parameter is a rough analogy and should be
refined against quantum-chemical calculations before quantitative use. Record the maximum
penalty in this table.

## The coordinates have to come too

The `itp` CHARMM-GUI writes fixes the order and the names of the atoms. **The
`<name>.mol2` and `<name>.pdb` that Open Babel made are in a different order and
must not be used for the simulation.** Put the coordinate file from the same
CHARMM-GUI job, usually `ligandrm.pdb`, in as `library/<name>/<RESI>.pdb`, or
`scripts/2_load.sh` stops before it starts.

Four files per molecule, then:

| file | where it comes from |
|---|---|
| `<RESI>.itp` | CHARMM-GUI Ligand Reader & Modeler, GROMACS output |
| `<RESI>.pdb` | the coordinates from that same job |
| `ff_<RESI>.itp` | the parameters that job supplied which the shipped force field lacks |
| `RESNAME.txt` | one line, `<RESI>` |

`ff_<RESI>.itp` is built from the job's `charmm36.itp`, keeping only the type
combinations `data/toppar/forcefield.itp` does not already define. Copying the
whole file instead redefines terms the polymer owns, and grompp refuses a second
block for the same dihedral type.

The `.rtf` and `.prm` from the job are worth keeping as well: they carry the
CGenFF penalties, and `f127load/check_params` reads them.
