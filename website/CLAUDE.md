# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 1. What this repository is

The `website/` folder of the GRABNES code repository (`github.com/gjung-group/grabnes`; formerly the separate repos `grabnes-docs` and `grabnes.github.com`) is the **source of the GRABNES project website**, published via GitHub Pages. The Fortran code is in the parent directory (`..`).

GRABNES (**GRA**phene/**B**N **E**lectronic **S**tructure) is a Fortran tight-binding package for layered 2D materials — graphene, hBN, twisted bilayer and multilayer graphene, graphene on hBN, etc. Key capabilities (from the old docs, `index.rst`):

- Builds commensurate moiré supercells with controllable twist angle and layer sliding (Hermann `n1 n2 n3 n4` indices); interfaces with LAMMPS for relaxation (`.mol` in / relaxed `.xyz` out).
- Observables: band structures, DOS (Lanczos recursion + exact diagonalization), unfolded spectral functions (ARPES comparison), LDOS, probability density, Kubo–Greenwood conductivity, magnetic-field (Landau-level) DOS with disorder.
- TB models: F2G2 (intralayer, from MLWFs), reparametrized two-center (TC) models for interlayer, "hybrid" F2G2+TC models, and the SHE model (TC variant tuned so the magic angle lands at a chosen value).
- OMP parallelization; tested with Intel compilers; built with a `make.sys` makefile (`make clean && make`, needs an empty `.version` file).

## 2. Commands

```bash
uv sync                                                   # create .venv from pyproject/uv.lock (MyST-MD needs Node.js >= 22)
uv run myst start                                         # live preview at http://localhost:3000/
uv run myst build --html --strict                         # what CI runs (static site in _build/html); must give zero warnings
uv run python tools/sync_from_code_repo.py [--check]      # check derived pages, regenerate API reference
uv run python tools/gen_api_reference.py [--src PATH] [--out PATH]   # regenerate reference/units.md + reference/api/ (sync runs it)
```

**The code repo is used in place, never copied.** This site is the `website/` folder of the code repository. The pages of `../docs/development/` and `../docs/user-guide/` are listed directly in `project.toc` (`file: ../docs/...`; MyST-MD gives pages outside the project folder flat URLs such as `/coding-conventions`), example figures and inputs are linked from `../examples/<ex>/`, and the old tutorial notebooks, figures and sample inputs from `../lanczosKuboCode/docs/{notebooks,figures}/`. Do not copy any of them into `website/`. Edit those files in the code repository; the theme's edit button points there (from `project.github` and the git path). `sync_from_code_repo.py` reads the HEAD commit of this checkout (`--src` for another): it regenerates the API reference (`reference/units.md`, `reference/api/`) when `lanczosKuboCode/Src` changed, and reports `MANUAL` for the hand-written pages derived from code-repo files (`development/index.md`, `theory/tapw-procedure.md`, `examples/*.md`, `user-guide/yaml-json-input.md`, and `installation.md`/`status.md`/`citing.md` from `README.md`, `lanczosKuboCode/README.md` and `CITATION.cff`), tracked by object hashes in `tools/sync-state.json`; update those pages by hand, verbatim, then re-run to record the new baseline. Last recorded: `1e390c6`.

There are no tests or linters; "passes" means `uv run myst build --html --strict` succeeds (the "Removing matplotlib string from outputs" notices are informational and do not fail it) and `uv run python tools/notebook_pages.py --check` reports no difference.

CI (`.github/workflows/docs.yml` at the code-repository root, steps run in `website/`, only when `website/**`, `docs/**`, `examples/**`, `lanczosKuboCode/docs/**` or the workflow changes): push to `main` → setup-node (24) → `uv sync` → `notebook_pages.py --check` → `BASE_URL=/<repository name> uv run myst build --html --strict --ci` → `upload-pages-artifact` (`website/_build/html`) → GitHub Pages (Settings → Pages → "GitHub Actions"). Also builds on PRs; deploys on push to `main` or `main-webpage` (the published site is built from `main-webpage`; the `github-pages` environment must allow that branch) or `workflow_dispatch`. Live site: <https://gjung-group.github.io/grabnes/>.

## 3. Current contents of this repo

```
README.md                  # intent of the site
grabnes.tar.gz             # 79 MB archive of the OLD sphinx tree — duplicate of grabnes/**, do not use
grabnes/var/www/sphinx/grabnesDocs/    # the OLD (2022) Sphinx docs — raw material ONLY
├── index.rst              # description, install, illustration gallery (some sections are empty stubs)
├── grabnes.rst            # Fortran API stubs (f:automodule — sources NOT in this repo)
├── tbModels.rst           # hybrid F2G2–TC model slides (images only, incomplete)
├── tutorial.rst, workshop.rst   # thin wrappers pointing at notebooks
├── conf.py                # old config (sphinxfortran + nbsphinx, classic theme)
├── figures/               # 10 PNGs
├── notebooks/             # the real content, see below
└── _build/                # stale generated HTML — never edit or reference
```

**Do not** build the old Sphinx tree in place, and do not add new pages to it. It is gitignored, read-only input for the migration. Never touch `grabnes.tar.gz` or `_build/` output; they will be removed once migration is done.

### The notebook content (primary source of prose)

- `notebooks/Twistronics2023.ipynb` — the flagship tutorial (~25k chars of markdown): graphene (build cell, bands, DOS) → AB-stacked bilayer & bulk graphite → twisted bilayer graphene (commensurate cells, Hamiltonian, Lanczos DOS, flat bands) → graphene on hBN (input flags, Hamiltonian visualization, magnetic-field DOS with/without Anderson disorder). Colab-oriented (`drive.mount`, Intel oneAPI setup); shown verbatim on the tutorial pages, with a warning on `tutorials/index.md` that it is a 2023 record.
- `notebooks/grabnes2.ipynb`, `grabnes3.ipynb` — input-file tutorials (`Gendata_graphene.in`, `Gendata_moire.in`, `ReadXYZ` for relaxed structures).
- `notebooks/results.ipynb` — plotting recipes (bands, DOS) in Python.
- Supporting files worth migrating next to the tutorial pages: `Gendata.in`, `BLBL.mol`, `BLBL.frac/.cart`, `moire.fdf`, `sublattices.dat`, `GenMoire.f90_tBG_template.f90`, `PyToGrabnes.py`, `createXYZ.py`, `fracToCart.py`, `replaceFrac.py`, `graphene/` sample outputs. Skip binaries/`a.out`, `__pycache__`, `.ipynb_checkpoints`, `.DS_Store`.

Known external links: Colab notebook `https://colab.research.google.com/drive/12_ofrjJi4pr6U9j-ER1WocYJkqxEffs-`, related code `https://github.com/gjung-group/JeilJungGroupCodes.git` (private: not linked on the site), `https://github.com/qtm-iisc/Twister`, ARPES paper `https://www.nature.com/articles/s41567-020-0974-x`. Author/copyright in old conf.py: Nicolas Leconte (2022), group: gjung-group (G. Jung, U. Liège).

## 4. Target website (agreed design decisions)

- **Stack: MyST-MD (Jupyter Book 2)** — the `mystmd` package (installed with uv; runs on Node.js), configured in `myst.yml` (project, table of contents, `book-theme` options). MyST Markdown pages and Jupyter notebooks (rendered natively, never executed: `jupyter: false`), KaTeX math, the `book-theme` layout with one stylesheet `assets/static/grabnes.css` (`site.options.style`) that only overrides the theme's `--myst-color-*` tokens with the logo palette (maroon/green/blue), adds gradient/hover accents (title, top bar, buttons, cards, active sidebar entry) and the page transition (a CSS entrance animation on `main.article-grid > article`, which the theme re-renders on every in-site navigation, plus `@view-transition` for full page loads; both off under `prefers-reduced-motion`). No custom JS or templates. The top bar (`site.nav`) has Overview (`/.`, the landing page `index.md`: MyST-MD refuses `/`, and `/index` is no client-side route, so clicking it shows "No result found for routeID"), Guide (`/guide`, `guide.md`: the documentation start page with the section cards) and Citing GRABNES (`/citing`; `citing.md` is `hidden` in the toc). Build via `myst build --html`; deploy with GitHub Actions to GitHub Pages. Moved from Sphinx/Furo in November 2026; page URLs have no `.html` suffix.
- **Notebooks → tutorial pages**: the notebooks are shown verbatim (not converted to prose), outputs included. The original `.ipynb` files are `../lanczosKuboCode/docs/notebooks/*.ipynb` of the code repository (downloads, not rendered), except `Twistronics2023.ipynb`, whose version with saved outputs lives in `tutorials/notebooks/` (a notebook found there is used instead of the code repo's). Every tutorial page `tutorials/<page>.ipynb` is a contiguous run of cells of one original notebook, copied by `tools/notebook_pages.py` (split defined in its `PAGES`; `--check` runs in CI) and rendered by MyST-MD with the saved outputs. The tool adds only a first front-matter cell (title, subtitle naming the notebook and its cells, `downloads` with the original `.ipynb` and the Colab link, `edit_url: null`) and two render fixes in the Markdown cells that change no text or formula (`fix_markdown`: the `PyToGrabnes.py` link target, and a blank line around a display `$$` block written directly under/above a text line, which MyST-MD otherwise mis-parses). Never hand-edit or paraphrase these pages, and do not restyle the cells.
- **Layout**: the **`website/` folder is the MyST project directory** — no `docs/` folder, no symlinks; every section is a real top-level folder:

```
myst.yml                   # MyST-MD config: project, toc (API part generated), book-theme options
index.md                   # Overview (landing page): what GRABNES is, capabilities, models, illustrations
guide.md                   # Guide: documentation start page (section cards, Colab getting started)
installation.md            # requirements, build, run, tests (quoted from the code repo README) + Colab
status.md                  # capabilities, support levels, known limitations (quoted from the README)
citing.md                  # citation, publications, authors, license (README + CITATION.cff)
pyproject.toml             # deps (see rule 7); uv sync
assets/                    # site-only files (code-repo files are linked in place, never copied)
├── images/                # site-only PNGs (model slides, code-structure diagram)
└── static/                # logo and favicon
theory/                    # models.md, moire-cells.md, observables.md, tapw-procedure.md
examples/                  # index + getting-started (01–04), moire (05–09), arpes (10–11, spectral functions)
user-guide/                # index, yaml-json-input (the other pages: ../docs/user-guide/)
tutorials/                 # introduction, graphene, bilayer-graphite, tbg,
                           # graphene-on-hbn, input-files, tbg-input, plotting
                           # (.ipynb cut from notebooks/ by tools/notebook_pages.py)
                           # + index; notebooks/ = Twistronics2023.ipynb (the others: ../lanczosKuboCode/docs/notebooks/)
benchmarks/                # benchmarks section sources
development/               # code-structure guide (index); the documents are ../docs/development/ of the code
                           # repo (functionality-status, validation, conventions…)
reference/index.md         # "API Reference" entry page (unit tables, links into units)
reference/units.md         # generated Fortran unit reference (overview) — never hand-edit
reference/api/             # generated: one page per Fortran unit, and api/<unit>/ one per public procedure — never hand-edit
tools/gen_api_reference.py # regenerates reference/units.md + reference/api/ (Ford + Doxygen comments)
```

  Non-content markdown (`README.md`, `CLAUDE.md`), and the old `grabnes/` tree are excluded in `myst.yml` (`project.exclude`). The navigation is `project.toc` in `myst.yml`: hand-edit it for new pages (no title-only groups: the theme carries their label into the previous/next links of every later page), except the block between the `BEGIN/END generated` markers, which `gen_api_reference.py` rewrites.

- Markdown pages cross-link with relative Markdown links (`[](../status.md)`, empty text = the target's title; `.ipynb` for tutorial pages) or `{ref}` to labels; buttons are the `{button}` role, cards `grid`/`grid-item-card` with `:link: page.md`, landing-page sub-lists the `{toc}` directive with `:context: children`; files of the code repository (figures, inputs, notebooks) are linked in place with relative `../` paths (see rule 3).
- Site tone: scientific documentation. Each tutorial page's download menu has "Open in Colab" and "Download <notebook>.ipynb" entries for its original notebook (front matter written by `notebook_pages.py`).
- "Edit this page" (the theme's edit button, from `project.github` in `myst.yml`): hand-written pages edit in this repo; pages of the code repo (`../docs/development/*`, `../docs/user-guide/*`) edit that file directly; generated pages (`reference/units.md`, `reference/api/`, `user-guide/input-keys.md`) and the tutorial notebook pages have no edit button (`edit_url: null`).

## 5. Migration map (old → new)

| Old source | New location | Notes |
|---|---|---|
| `index.rst` "Description" | `index.md` | rewrite, add feature cards (grid) |
| `index.rst` "Getting started"/"Installation" | `installation.md` | keep Colab + compile paths |
| `index.rst` "Illustrations" (tBG, t20G, ARPES) | `benchmarks/index.md` (+ pages) | empty stubs (gap estimation, supermoire, impurities) → marked *coming soon*, do not invent content |
| `tbModels.rst` + `GBNtwoLayers*.png` | `theory/models.md` | slides are image-only; write prose around them from notebook theory |
| `grabnes.rst` (f:automodule stubs) | `reference/index.md` | module list + `codeStructure.png` + link to code repo; sphinxfortran is dead, drop it |
| `Twistronics2023.ipynb` §Introduction + prerequisites | `tutorials/introduction.ipynb` | verbatim cells (`tools/notebook_pages.py`) |
| `Twistronics2023.ipynb` §Graphene | `tutorials/graphene.ipynb` | verbatim cells |
| `Twistronics2023.ipynb` §Bernal/graphite | `tutorials/bilayer-graphite.ipynb` | verbatim cells |
| `Twistronics2023.ipynb` §tBG | `tutorials/tbg.ipynb` | verbatim cells |
| `Twistronics2023.ipynb` §Graphene on hBN + magnetic DOS + Future work | `tutorials/graphene-on-hbn.ipynb` | verbatim cells |
| `grabnes2.ipynb` | `tutorials/input-files.ipynb` | whole notebook, verbatim |
| `grabnes3.ipynb` | `tutorials/tbg-input.ipynb` | whole notebook, verbatim |
| `results.ipynb` | `tutorials/plotting.ipynb` | whole notebook, verbatim |
| `figures/*.png` | `../lanczosKuboCode/docs/figures/` in place; the others `assets/images/` | |
| `notebooks/*.{in,mol,dat,f90,py}` | `../lanczosKuboCode/docs/notebooks/` in place | sample inputs/scripts |
| code repo `examples/*/README.md` + figures | `examples/{getting-started,moire,arpes}.md` | figures and inputs linked in place from `../examples/` |
| code repo `docs/user-guide/*.md` | in place (toc) | `input-keys.md` is generated from source — regenerate, never hand-edit |
| code repo `docs/development/*.md` | in place (toc) | `README.md` folded into `development/index.md` |
| code repo `docs/theory/tapw-procedure.tex` | `theory/tapw-procedure.md` | LaTeX → MyST, formulas verbatim |

## 6. Working rules

1. **Never fabricate or modify existing information for any system.** All theory, parameters, formulas, flags, and code snippets must appear exactly as in the old docs/notebooks (reformatting to MyST is fine; changing values, formulas, or comments is not). If a source looks like a typo (e.g. the `i^2+i^2+ij` denominator in the α formula), keep the original verbatim and flag it in an admonition — never silently "fix" it. Empty stub sections stay marked "coming soon".
2. Markdown (MyST) for all new pages; LaTeX math via `$...$`/`$$...$$`; cross-page links via relative `.md` links or `{ref}`; headings must keep implicit slugs stable once published.
3. Images: a figure of the code repository is linked in place (`../examples/<ex>/reference/...`, `../lanczosKuboCode/docs/figures/...`, one more `../` from section pages); only site-only images go in `assets/images/`. Always with alt text. No base64 blobs in Markdown.
4. Keep pages buildable: `uv run myst build --html --strict` (from the repo root) must succeed with zero warnings before committing.
5. Do not commit `_build/`, `__pycache__`, `.ipynb_checkpoints`, `.DS_Store` (see `.gitignore`).
6. The old tree under `grabnes/` stays read-only until migration is verified; delete it (and `grabnes.tar.gz`) only in a final, separate commit.
7. Python tooling: declare deps in `pyproject.toml` (`mystmd`, `ford`, `gitpython`); `uv sync` builds `.venv` from `uv.lock`, and CI installs the same way (plus Node.js 24 via `actions/setup-node`).
8. `reference/units.md` and `reference/api/**/*.md` are generated (Ford over the HEAD of this checkout, run on a comment-normalized temp copy, plus the Doxygen blocks (`@brief`, `@param`, `@par`, `@see`, ...) or, without one, the plain comment above each definition (MIO `INPUT`/`OUTPUT` headers give the parameter descriptions), and the `!! - \`Key\`:` input-key items of the sources, all quoted verbatim; calls/callers from Ford's `correlate()`): after upstream changes run `uv run python tools/sync_from_code_repo.py` and commit the result; never hand-edit these files.
9. Commits: small, imperative-mood messages ("Add theory/models page…").

## 7. Open questions / follow-ups

- [ ] Canonical URL of the Fortran code repo (for Installation/Reference links) — likely `github.com/gjung-group/grabnes`; confirm before linking.
- [x] `README.md` mentioned `grabnes/examples/*` — resolved: the example pages `examples/*.md` link the inputs and figures of the code repo's `examples/` in place.
- [ ] Fill empty illustration stubs (G/hBN gap, supermoire, impurity states) — needs content from maintainers.
- [x] Logo/branding: `assets/static/logo1.png` (wordmark banner, added by the maintainers) and `logo-icon.png` (square icon cropped from it, also the site favicon); referenced from `README.md` and set as `site.options.logo`/`favicon` in `myst.yml`.
- [ ] Remove `grabnes/` and `grabnes.tar.gz` after final migration commit.
