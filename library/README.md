# Guest library

Molecules that are commonly loaded into Pluronic F127 micelles, prepared for use with this
workflow. Properties are taken from PubChem (CID given); XLogP is the PubChem computed value.

3D structures (`*.mol2`, `*.pdb`) were generated from the canonical SMILES with Open Babel
(`--gen3d --minimize --ff MMFF94`) and are provided **only as input for parameter
generation**. They are not equilibrated.

## Status

| Molecule | CID | Formula | MW | XLogP | Atoms | Topology | Tier |
|---|---|---|---|---|---|---|---|
| pyrene       | 31423   | C16H10     | 202.3 | 4.9 | 26  | ✅ `PYR.itp` | A |
| ibuprofen    | 3672    | C13H18O2   | 206.3 | 3.5 | 33  | ⬜ | A |
| aspirin      | 2244    | C9H8O4     | 180.2 | 1.2 | 21  | ⬜ | A |
| curcumin     | 969516  | C21H20O6   | 368.4 | 3.2 | 47  | ⬜ | A |
| indomethacin | 3715    | C19H16ClNO4| 357.8 | 4.3 | 41  | ⬜ | A |
| nile_red     | 65182   | C20H18N2O2 | 318.4 | 3.8 | 42  | ⬜ | A |
| resveratrol  | 445154  | C14H12O3   | 228.2 | 3.1 | 29  | ⬜ | A |
| quercetin    | 5280343 | C15H10O7   | 302.2 | 1.5 | 32  | ⬜ | B |
| paclitaxel   | 36314   | C47H51NO14 | 853.9 | 2.5 | 113 | ⬜ | B |
| docetaxel    | 148124  | C43H53NO14 | 807.9 | 1.6 | 111 | ⬜ | B |
| doxorubicin  | 31703   | C27H29NO11 | 543.5 | 1.3 | 68  | ⬜ | C |

**Tier A** — small, neutral, drug-like. CGenFF penalties are expected to be low.
**Tier B** — large or highly substituted. Check the penalty report; some dihedrals may need
attention.
**Tier C** — **not supported as-is.** Doxorubicin carries a protonated amine at physiological
pH and is therefore cationic. The loading mechanism described in the accompanying paper is
purely dispersive, so a charged guest is expected to behave differently and is included here
only as a documented limitation.

Note that the carboxylic acids (ibuprofen, aspirin, indomethacin) are supplied as the
**neutral** forms. At pH 7.4 they are largely deprotonated; the neutral form is the species
that partitions into the core, but this should be stated explicitly in any work using them.

## Adding a topology

1. Open [CHARMM-GUI Ligand Reader & Modeler](https://charmm-gui.org/?doc=input/ligandrm).
2. Upload `library/<name>/<name>.mol2` (or paste the SMILES from `SMILES.txt`).
3. Select **GROMACS** output. Download the job archive.
4. From the archive, copy into `library/<name>/`:
   * `toppar/<RESI>.itp` → the guest topology
   * the CGenFF penalty report, saved as `penalty.txt`
5. Create `RESNAME.txt` containing the residue name used in the itp (e.g. `IBU`).
6. Record the CHARMM-GUI job id in `CHARMMGUI_JOBID.txt` for reproducibility.
7. Verify with `bash tests/smoke_test.sh <name>`.

**Check the penalty report.** CGenFF reports a penalty for every assigned parameter. Values
below 10 are generally acceptable; above 50 the parameter is a rough analogy and should be
refined against quantum-chemical calculations before quantitative use. Record the maximum
penalty in this table.

## 좌표 파일도 함께 넣어야 한다

CHARMM-GUI 가 만든 `itp` 는 원자 순서와 이름이 정해져 있다. **Open Babel 이 만든
`<name>.mol2` / `<name>.pdb` 는 그 순서와 다르므로 시뮬레이션에 쓰면 안 된다.**
CHARMM-GUI 작업 결과에 들어 있는 좌표 파일(보통 `<RESI>.pdb` 또는 `<RESI>.crd` → pdb 변환본)을
`library/<name>/<RESI>.pdb` 로 함께 넣어야 `scripts/2_load.sh` 가 동작한다.

정리하면 분자 하나당 필요한 파일은 셋이다.

| 파일 | 출처 |
|---|---|
| `<RESI>.itp` | CHARMM-GUI Ligand Reader & Modeler (GROMACS 출력) |
| `<RESI>.pdb` | 같은 작업의 좌표 파일 |
| `RESNAME.txt` | `<RESI>` 한 줄 |

`penalty.txt` 와 `CHARMMGUI_JOBID.txt` 는 재현성을 위해 권장한다.
