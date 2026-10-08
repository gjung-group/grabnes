# Third-party components and notices

GRABNES as a whole is distributed under the GNU General Public License,
version 3 or (at your option) any later version; see [`LICENSE`](LICENSE).
The components below keep the copyright and license notices of their authors.
Those notices are not replaced by the project license: each file remains
under the terms stated in it, and the combination is distributed under the
GPL because all of those terms permit it.

| Component | Files | Notice in the files | Terms |
| --- | --- | --- | --- |
| MIO library (memory, input/output, timers, MPI wrappers) | `lanczosKuboCode/Src/MIO/`, `lanczosKuboCode/Src/MIO/MPI/` | Copyright (C) 2012 Rafael Martinez-Gordillo; "distributed under the terms of the GNU General Public License" | GNU GPL, no version specified. Under the license's own terms a recipient may then choose any published version; it is used here under version 3 or later. License text: `LICENSE` |
| Math library | `lanczosKuboCode/Src/math/` | same notice | same |
| Solver build file | `lanczosKuboCode/Src/Makefile` | "Kubo3 Makefile, Copyright (C) 2012 Rafael Martinez-Gordillo" | no separate terms stated; distributed with the solver |
| `cosine_transform_data`, `cosine_transform_inverse`, `r8vec_uniform_01`, `timestamp` | end of `lanczosKuboCode/Src/ham.F90` | Author: John Burkardt; "distributed under the GNU LGPL license" | GNU Lesser GPL, no version specified. The LGPL permits use in a GPL-licensed program. License text: `COPYING.LESSER` (version 3), which supplements `LICENSE` |
| `fracToCart.py` | eleven copies under `LAMMPSNotebooks/` and `PyBinding/` | Copyright (c) 2014, Gavin Heverly-Coulson, with the conditions and disclaimer of a two-clause BSD license in each file | Two-clause BSD license as printed in the files; the notice must stay with them |

Libraries that are linked at build time but not distributed with GRABNES:

| Library | Use | License of the library |
| --- | --- | --- |
| BLAS, LAPACK | dense linear algebra | modified BSD (reference implementation); other implementations have their own terms |
| ARPACK | sparse eigenvalue problems | BSD-style |
| An MPI implementation | program start-up and communication layer | depends on the implementation |

A binary that includes such a library must also satisfy that library's
license.

Removed: two commented-out routines taken from *Numerical Recipes* (`tqli`,
`pythag`) were deleted from `lanczosKuboCode/Src/kubo.F90`. They were never
compiled, and their source may not be redistributed.
