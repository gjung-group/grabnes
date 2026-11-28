# GRABNES website

<p align="center">
  <img src="assets/static/logo1.png" alt="GRABNES logo" width="480">
</p>

Source of the documentation website for **GRABNES**, published with GitHub
Pages at <https://gjung-group.github.io/grabnes/>. It lives in the `website/`
folder of the [code repository](https://github.com/gjung-group/grabnes).

## About GRABNES

GRABNES (**GRA**phene/**B**N **E**lectronic **S**tructure) is a Fortran
tight-binding package from the
[gjung-group](https://github.com/gjung-group) for layered 2D materials:
graphene, hBN, twisted bilayer and multilayer graphene, graphene on hBN, and
related heterostructures.

- **Moiré supercells:** builds commensurate cells with controllable twist
  angle and layer sliding, and interfaces with LAMMPS for structural
  relaxation.
- **Observables:** band structures, density of states (Lanczos recursion and
  exact diagonalization), unfolded spectral functions for comparison with
  ARPES, LDOS, probability densities, Kubo–Greenwood conductivity, and
  magnetic-field (Landau-level) DOS with disorder.
- **Tight-binding models:** F2G2 intralayer models from maximally localized
  Wannier functions, reparametrized two-center (TC) interlayer models, hybrid
  F2G2+TC models, and the SHE model.
- **Performance:** OpenMP parallelization.

## What the website covers

- Installation: compiling the code, or running it in the browser via Colab.
- Theory of the implemented models and observables.
- Worked examples (the supported examples 01–11 of the code repository) and
  step-by-step tutorials.
- A user guide: the input file, writing it in YAML or JSON, stacking
  displacements, TAPW + SOC input, and the complete input-key reference with
  a description of every key.
- Benchmark studies (some still marked "coming soon").
- A development guide and a generated Fortran API reference: one page per
  unit with signatures, arguments, descriptions from the source comments,
  dependencies, source links, and the input keys each routine reads.

## Current state

The site is a [MyST-MD](https://mystmd.org) (Jupyter Book 2) project whose
**source directory is `website/`**: `myst.yml` holds the configuration
and the table of contents, the pages are MyST Markdown and Jupyter notebooks
(rendered natively with their saved outputs, not executed), and the look is
MyST's `book-theme` as it ships (single-page navigation, hover previews of
cross-references, search, light and dark mode). It is deployed to GitHub Pages
by a GitHub Actions workflow.
Content comes from two sources: the old (2022) Sphinx docs and notebooks
(kept locally under the gitignored `grabnes/` tree as read-only migration
input, to be removed once the migration is verified), and the files of the
[GRABNES code repository](https://github.com/gjung-group/grabnes) that contains
this folder, used in place: `docs/development/`, `docs/user-guide/`, `examples/`
and `lanczosKuboCode/docs/` (notebooks, figures, sample inputs) are listed in
the table of contents or linked with `../` paths, never copied here.

Sections:

| Folder | Content |
| --- | --- |
| `myst.yml` | MyST-MD configuration: project, table of contents (the API part is written by `gen_api_reference.py`), theme options |
| `index.md`, `installation.md` | landing page; requirements, build, run and tests (quoted from the code repository README) and the 2023 Colab notebook |
| `status.md`, `citing.md` | supported functionality and known limitations; citation, authors and license (quoted from the README and `CITATION.cff`) |
| `theory/` | tight-binding models, moiré supercells, observables, TAPW procedure |
| `examples/` | the supported examples of the code repository (01–11: getting started, moiré, ARPES) with figures and downloadable inputs |
| `tutorials/` | the cells of the original Jupyter notebooks, shown unchanged with their saved outputs (not executed, rendered by MyST-MD): Twistronics 2023 workshop (introduction, graphene, bilayer/graphite, twisted bilayer, graphene on hBN) and the moiré-course notebooks (input files, tBG input files, plotting); the originals are `lanczosKuboCode/docs/notebooks/*.ipynb` of the code repository, except `Twistronics2023.ipynb`, whose version with saved outputs is kept in `notebooks/` |
| `user-guide/` | index and the YAML/JSON input converter page; the other user-guide pages are `docs/user-guide/*.md` of the code repository |
| `benchmarks/` | benchmark studies; empty illustration stubs are marked "coming soon" |
| `development/` | index (code structure); the other development pages are `docs/development/**.md` of the code repository |
| `reference/` | API reference entry page + generated Fortran unit reference (overview `units.md`, one page per unit in `api/`) |
| `assets/` | images that exist only on the site, logo |
| `tools/` | check of the pages derived from the code repository (`sync_from_code_repo.py`, state in `sync-state.json`) and the API-reference generator (`gen_api_reference.py`), and `notebook_pages.py`, which cuts the tutorial pages out of the original notebooks |

## Keeping in step with the code repository

Nothing is mirrored: pages, figures and inputs of the code repository are used
in place, so they are always current. `reference/units.md` and `reference/api/` are generated from its sources with
Ford (one page per unit and one per public procedure), plus the Doxygen
(`@brief`, `@param`, ...) blocks or plain comments above each definition and the
`!!` input-key comments of the sources, quoted verbatim. Never hand-edit the generated reference.

```bash
uv run python tools/sync_from_code_repo.py [--check]   # HEAD of this checkout; --check only reports drift (exit 1)
uv run python tools/gen_api_reference.py               # regenerate reference/units.md + reference/api/
```

Hand-written pages derived from upstream are tracked in
`tools/sync-state.json`; the script prints `MANUAL` when their source changed.

The script is run by hand; there is no automatic workflow. Commit the
result and push to `main` to deploy it.

## Building locally

```bash
uv sync                                        # creates .venv from pyproject/uv.lock (needs Node.js >= 22)
uv run myst start                              # live preview at http://localhost:3000/
uv run myst build --html                       # static site in _build/html/
```

The tutorial pages `tutorials/*.ipynb` are copied from the original notebooks
(`../lanczosKuboCode/docs/notebooks/`, or `tutorials/notebooks/` for the site's own version), with a generated front matter cell in front (title,
source cells, Colab and download links) and whitespace-only render fixes
(listed in the tool); never edit them by hand. After changing a notebook
or the page split (`PAGES` in `tools/notebook_pages.py`), run
`uv run python tools/notebook_pages.py`; CI runs it with `--check` and fails
if a page differs from its notebook cells.

The CI build treats warnings as errors; check with
`uv run myst build --html --strict` before committing.

## Deployment

The live site is at <https://gjung-group.github.io/grabnes/>.
`.github/workflows/docs.yml` at the root of the code repository builds the site to `website/_build/html/`
(with `BASE_URL=/<repository name>`) on every push to
`main-webpage` (the branch the published site is built from) or `main`, or when run by hand
from the Actions tab, and deploys it to GitHub Pages (Settings → Pages → "GitHub Actions";
the `github-pages` environment must allow the branch). Pull requests are built too
(with `--strict`) but not deployed. Page URLs have no `.html` suffix
(e.g. `reference/api/sys/syskill`).
