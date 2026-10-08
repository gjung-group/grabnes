# Licensing and third-party code audit

Audit of the tracked files as of October 2026, and the license decision taken
on its basis. It records what the files say and is not legal advice.

## Decision

GRABNES is distributed under the **GNU General Public License, version 3 or
any later version** (`GPL-3.0-or-later`). The maintainer selected it on
2026-10-09. The license text is in `LICENSE`; the notices of other authors
are listed in `THIRD_PARTY_LICENSES.md` and remain in their files.

Why this license is compatible with what the repository contains:

1. **GPL notices without a version (MIO and math libraries, 33 files).**
   The notices say "distributed under the terms of the GNU General Public
   License" and name no version. The license itself settles this case: "If
   the Program does not specify a version number of the GNU General Public
   License, you may choose any version ever published by the Free Software
   Foundation" (GPL version 3, section 14; versions 1 and 2 contain the same
   rule in sections 7 and 9). The libraries may therefore be used and
   redistributed under version 3 or later. Their notices are unchanged; the
   `LICENSE` file they refer to now exists.
2. **LGPL routines (four routines by John Burkardt in `ham.F90`).** The LGPL
   allows the covered code to be combined with, and conveyed as part of, a
   GPL-licensed program. The routines keep their LGPL notice; the LGPL text
   is supplied as `COPYING.LESSER`.
3. **Two-clause BSD script (`fracToCart.py`).** This license is compatible
   with the GPL; its notice and disclaimer stay in each copy.
4. **Solver sources without a notice.** They were written by the three
   authors of GRABNES (Rafael Martinez-Gordillo's original program, extended
   in the group of Jeil Jung, currently developed by Nicolas Leconte), with
   contributions from other group members named in the sources. The original
   program is the one whose libraries carry the GPL notice, so a GPL license
   for the whole is consistent with the terms under which it was received.
5. **A copyleft license is the conservative choice.** A permissive license
   for the whole would have required the copyright holder of the GPL-marked
   libraries to relicense them; the GPL does not.

What was **not** done: no copyright notice was removed or reworded, no file
was relicensed away from the terms stated in it, and no copyright ownership
was asserted on behalf of anyone. Authorship (who wrote the software) is
recorded in `CITATION.cff`; it is not a statement about who holds copyright.

## Open points

| Point | Whose confirmation | Files |
| --- | --- | --- |
| Agreement of the co-authors with `GPL-3.0-or-later` for the solver sources that carry no notice | Jeil Jung, Rafael Martinez-Gordillo | `lanczosKuboCode/Src/*.F90`, `*.f90`, `Utils/`, the Makefiles |
| Whether copyright headers are added to those sources, and naming whom | the three authors and, where applicable, their institutions | same |
| Terms for the research material | maintainer | `LAMMPSNotebooks/`, `PyBinding/`, `usefulGeneralScripts/`: notebooks, data and scripts without notices, currently covered by the repository license by default |
| Contributions of other group members named in source comments | maintainer | `ham.F90` and others |

None of these prevents distribution under the GPL as it stands; they are the
points on which a written record would be prudent before a formal release.

## Findings of the audit

| Component | Location | Notice found | Remark |
| --- | --- | --- | --- |
| Project metadata (before this change) | `pyproject.toml` | `license = {text = "Proprietary"}` | Contradicted a public repository and the notices below |
| Top level (before this change) | | no `LICENSE` or `COPYING` file | The GPL notices below refer to "the file `LICENSE` in the root directory" |
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

## Changes made after the audit

- `LICENSE` (GPL version 3) and `COPYING.LESSER` (LGPL version 3) added.
- `THIRD_PARTY_LICENSES.md` added.
- `pyproject.toml`: "Proprietary" replaced by `GPL-3.0-or-later`; authors and
  maintainer recorded.
- The two commented-out *Numerical Recipes* routines (`pythag`, `tqli`; 77
  comment lines) were removed from `lanczosKuboCode/Src/kubo.F90`. They were
  comments only, referenced nowhere in compiled code, and their removal
  changes no executable statement. The group's own commented routine that
  called them was left in place.
