# Licensing and third-party code audit

Audit of the tracked files as of October 2026. It records what the files say;
it is not legal advice, and **no license was selected, added, or changed**.

## Findings

| Component | Location | Notice found | Remark |
| --- | --- | --- | --- |
| Project metadata | `pyproject.toml` | `license = {text = "Proprietary"}` | Contradicts a public repository and the notices below |
| Top level | | no `LICENSE` or `COPYING` file | The GPL notices below refer to "the file `LICENSE` in the root directory", which does not exist |
| MIO library (memory, input/output, timers, MPI wrappers) | `lanczosKuboCode/Src/MIO/`, `Src/MIO/MPI/` (25 files) | "Copyright (C) 2012 Rafael Martinez-Gordillo ... distributed under the terms of the GNU General Public License" | No GPL version is named; the notice points to `http://www.gnu.org/copyleft/gpl.txt` |
| Math library | `lanczosKuboCode/Src/math/` (8 files) | same notice | same |
| Build file of the solver | `lanczosKuboCode/Src/Makefile` | "Kubo3 Makefile, Copyright (C) 2012 Rafael Martinez-Gordillo" | No license statement |
| Solver sources | `lanczosKuboCode/Src/*.F90`, `*.f90` (27 files) | no copyright or license header | `ham.F90` names as authors "Rafael Martinez-Gordillo (original), Nicolas (modifications), Jinwoo (current version)". The solver descends from the same 2012 code ("Kubo3") as the two libraries |
| Four utility routines | end of `lanczosKuboCode/Src/ham.F90` (`cosine_transform_data`, `cosine_transform_inverse`, `r8vec_uniform_01`, `timestamp`) | "Author: John Burkardt ... distributed under the GNU LGPL license" | Third-party code; notices must be kept |
| Commented-out routines | `lanczosKuboCode/Src/kubo.F90`, about lines 420 to 540 | comment "Fait appel a tqli.f, donc aussi a pythag.f (Numerical Recipes)" followed by the routines as comments | Not compiled, but the text of Numerical Recipes routines is not freely redistributable |
| `fracToCart.py` (eleven copies) | `LAMMPSNotebooks/`, `PyBinding/` | "Copyright (c) 2014, Gavin Heverly-Coulson" with a two-clause BSD-style text | Third-party script; notice must be kept with the file |
| Sphinx configuration | `docs/conf.py`, `lanczosKuboCode/docs/conf.py` | author "Jeil Jung Group"; copyright "2022, Nicolas Leconte" | Metadata only |
| Notebooks, data, scripts | `LAMMPSNotebooks/`, `PyBinding/`, `usefulGeneralScripts/` | none | Authorship and reuse terms not stated |

Libraries that are linked but not distributed with GRABNES (BLAS, LAPACK,
ARPACK, MPI) impose no condition on the source distribution; their own
licenses apply to binaries that include them.

## What this implies

1. **The GPL-marked libraries are part of the executable.** MIO and the math
   library are compiled into `grabnes`. As long as they are under the GPL, the
   program as a whole can only be distributed under GPL-compatible terms,
   unless their copyright holder agrees to other terms. "Proprietary" in
   `pyproject.toml` is not compatible with distributing them.
2. **The GPL version is unspecified.** The notices name no version. This needs
   to be settled with the copyright holder (GPL-2.0-only, GPL-2.0-or-later, and
   GPL-3.0 differ in their compatibility with other licenses).
3. **The origin of the unmarked solver sources should be recorded.** They
   derive from the same original code; whoever holds rights in the 2012 code
   has to agree to the license chosen for the whole.
4. **LGPL and BSD-style notices must be preserved** in the files that carry
   them, whatever project license is chosen.
5. **The commented-out Numerical Recipes routines** are the one item that
   could not be redistributed under any open license. They are dead code.

## Decisions required from the maintainers

1. Choose the project license, in agreement with Rafael Martinez-Gordillo as
   the named copyright holder of MIO and the math library (and, as far as it
   applies, of the original solver).
2. State the GPL version of the existing notices, or have them replaced by
   the copyright holder.
3. Add the corresponding `LICENSE` file at the top level and make
   `pyproject.toml` agree with it.
4. Decide who is named as copyright holder of the solver sources, and whether
   headers are added to them.
5. Authorize the removal of the commented-out Numerical Recipes routines from
   `kubo.F90`. (Not done here: `kubo.F90` is a scientific source file under
   the current code freeze, and the removal is a licensing decision.)
6. Decide whether `LAMMPSNotebooks/`, `PyBinding/`, and
   `usefulGeneralScripts/` are released under the same terms, under separate
   terms, or moved out of the software distribution.

Nothing in this list was acted upon. No copyright notice was removed or
altered during the release preparation.
