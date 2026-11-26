# MoleKing backlog for the EField skill

Written 2026-11-26 from a read of `main` (0ece191). Branch for this work: `Revision`.
The EField skill will generate every Gaussian/ORCA input and read every output through MoleKing.
Items marked **(verify)** need a check against the program manual before coding.

## 1. Alignment (needed by all field protocols)

1. **Three-atom frame.** Add e.g. `alignFrame(a, b, c, axis='z', plane='xz')`: atom a at the origin,
   a→b along +axis, atom c in the chosen plane on its positive side. Deterministic, documented, returns
   (or stores) the rotation matrix applied. Today `alignBond(a, b)` leaves the rotation about the axis
   free, and the `toGJF(zmatrix=True, bondFixer=[a,b])` frame depends implicitly on the third atom of
   the reordered Z-matrix.
2. **Vector transform helpers.** Rotate a field vector between the input frame and the molecular frame,
   so a field can be specified in the molecular frame and written in the input frame (or the reverse).

## 2. `Molecule.toGJF`

3. **Bug:** `else if (bondFixer.size() == 1 || bondFixer.size() > 2)` is unreachable (the first branch
   catches any non-empty list), so invalid `bondFixer` sizes never raise.
4. **Field with Cartesian coordinates.** The warning "Geometry Optimization with Electric field can only
   be applied to zmatrix" forces Z-matrix inputs (also for TS searches). Allow Cartesian + `NoSymm`
   optimisations with `Field=Read` **(verify** that Gaussian keeps the input orientation**)**.
5. **`NoSymm` added automatically** whenever `EField` is set (check current behaviour).
6. **Field precision.** `to_string` writes 6 decimals; allow a format (e.g. `%.8e`) for small fields.
7. **Frozen atoms.** Option to freeze by index list or by element (`opt=ReadOpt` with
   `noatoms atoms=... notatoms=...`, or `-1` flags), and return the frozen list so it can be logged and
   checked.
8. **Per-element basis/ECP** (`Gen`/`GenECP` blocks, e.g. 6-311+G(d,p) on light atoms + SDD on Au).
9. **Multi-step jobs.** `--Link1--` steps (opt+freq → single point, `field=checkpoint`,
   `geom=check guess=read`, `%oldchk`).
10. **Keywords as arguments:** `Temperature=`, `Stable=Opt`, `Opt=(TS,CalcFC,NoEigenTest)`.

## 3. ORCA

11. **Merge the ORCA parser** (`origin/ORCAlogfile`) into `main`; `main` has no ORCA reader.
12. **`toORCA()` input writer:** method/basis, `%pal`, `%maxcore`, charge/multiplicity, Cartesian
    coordinates, per-element basis/ECP, frozen atoms (`%geom Constraints {C n C} end`), TS optimisation
    (`OptTS`, initial Hessian), `Freq`, stability analysis, point charges (external file), and the field
    (`%scf EField x, y, z end`) **(verify** syntax, units, sign convention, and whether gradients with
    a field are supported for optimisation in the installed ORCA version**)**.
13. **ORCA parser fields:** termination, final energy, dipole vector, frequencies (imaginary count),
    thermochemistry, atomic charges, printed field, final geometry, stability result.

## 4. `G16LOGfile` (fields the skill's checker needs)

14. Termination status per link (Normal / Error, error message).
15. Imaginary-frequency count (helper over `getVibFrequencies`).
16. Route section per link; printed field values (to compare with the planned field).
17. Choice of orientation in `getMolecule` (input vs standard).
18. Dipole vector per link, with its orientation stated.
19. Atomic charges (Mulliken, Hirshfeld) per atom, for group charges.
20. Wavefunction stability result.

## 5. Tests

21. Tests for the three-atom frame, field writing (Gaussian and ORCA), frozen-atom writing, and the new
    parser fields, using small molecules.

## 6. Repository note

22. `.git/sequencer` holds an unfinished multi-commit `git revert` (head 73bbd43; commits f3627d1 …
    0ed1746 on pyproject.toml/setup.py). `git status` reports "Revert currently in progress" and refuses
    to switch branches. Resolve with `git revert --quit` (keeps the current HEAD) or `git revert --abort`.
    The `Revision` branch was created with `git branch` (no checkout) for this reason.
