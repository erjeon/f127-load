# Solute library

One directory per molecule. Properties are from PubChem (CID given); XLogP is the PubChem
computed value. The `*.mol2` and `*.pdb` made from SMILES with Open Babel are input for
parameter generation only.

## Status

| molecule | CID | formula | MW | XLogP | atoms | charge | topology |
|---|---|---|---|---|---|---|---|
| pyrene            | 31423   | C16H10      | 202.3 | 4.9 | 26  |  0 | `PYR.itp` |
| aspirin           | 2244    | C9H8O4      | 180.2 | 1.2 | 21  |  0 | `ASP.itp` |
| curcumin          | 969516  | C21H20O6    | 368.4 | 3.2 | 47  |  0 | `CUR.itp` |
| paclitaxel        | 36314   | C47H51NO14  | 853.9 | 2.5 | 113 |  0 | `TAX.itp` |
| doxorubicin       | 31703   | C27H29NO11  | 543.5 | 1.3 | 68  |  0 | `DOX.itp` |
| doxorubicin (+1)  | 31703   | C27H30NO11+ | 544.5 | 1.3 | 69  | +1 | `DXP.itp` |
| ibuprofen         | 3672    | C13H18O2    | 206.3 | 3.5 | 33  |  0 | structure only |
| indomethacin      | 3715    | C19H16ClNO4 | 357.8 | 4.3 | 41  |  0 | structure only |
| nile_red          | 65182   | C20H18N2O2  | 318.4 | 3.8 | 42  |  0 | structure only |
| resveratrol       | 445154  | C14H12O3    | 228.2 | 3.1 | 29  |  0 | structure only |
| quercetin         | 5280343 | C15H10O7    | 302.2 | 1.5 | 32  |  0 | structure only |
| docetaxel         | 148124  | C43H53NO14  | 807.9 | 1.6 | 111 |  0 | structure only |

Only molecules with a topology are offered by `./f127 new` and the browser page.
`python3 scripts/check_library.py` prints this state for every entry.

Doxorubicin ships twice: `DOX` is the neutral molecule as CHARMM-GUI returns it, `DXP`
carries the protonated daunosamine amine (net charge +1), the species at physiological pH.
The carboxylic acids (ibuprofen, aspirin, indomethacin) are supplied as the neutral forms.

## Adding a topology

1. Open [CHARMM-GUI Ligand Reader & Modeler](https://charmm-gui.org/?doc=input/ligandrm).
2. Upload `library/<name>/<name>.mol2`, or paste the SMILES from `SMILES.txt`.
3. Select GROMACS output and download the job archive.
4. Copy into `library/<name>/`:
   * `toppar/<RESI>.itp` as `<RESI>.itp`
   * `ligandrm.pdb` as `<RESI>.pdb` (the coordinates from the same job; the Open Babel
     files are in a different atom order and are not used for the simulation)
   * the `.rtf` and `.prm` from the job, which carry the CGenFF penalties
5. Write `RESNAME.txt` with the three-letter residue name used in the itp.
6. Build the force field fragment and check the entry:

   ```bash
   python3 scripts/make_ff_fragment.py library/<name> /path/to/charmm36.ff
   python3 scripts/check_library.py <name>
   python3 -m f127load.check_params library/<name>
   bash tests/smoke_test.sh <name>
   ```

`ff_<RESI>.itp` keeps only the type combinations that `data/toppar/forcefield.itp` does not
already define. CGenFF penalties below 10 are fine; above 50 the parameter should be refined
before quantitative use.

| file | from |
|---|---|
| `<RESI>.itp` | CHARMM-GUI Ligand Reader & Modeler, GROMACS output |
| `<RESI>.pdb` | the coordinates from that same job |
| `ff_<RESI>.itp` | `make_ff_fragment.py` |
| `RESNAME.txt` | one line, `<RESI>` |
