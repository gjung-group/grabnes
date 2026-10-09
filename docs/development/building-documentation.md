# Building the GRABNES documentation

GRABNES combines Doxygen's Fortran parser with Sphinx and Breathe. Doxygen
extracts API information from `lanczosKuboCode/Src`; Sphinx renders that API
alongside the Markdown guides under `docs/`.

## Requirements

- Doxygen
- Python 3.11 or later
- the packages pinned in `docs/requirements.txt`

Use an isolated Python environment if possible:

```sh
python3 -m venv .venv-docs
. .venv-docs/bin/activate
python3 -m pip install -r docs/requirements.txt
```

The local environment and generated `_build` directory must not be committed.

## Local build

Run these commands from the repository root:

```sh
doxygen Doxyfile
python3 -m sphinx -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser. `-W` makes Sphinx report
documentation warnings as build failures while `--keep-going` displays all
problems in one run.

## Configuration map

- `Doxyfile` selects the maintained Fortran source tree and writes XML to
  `docs/_build/doxygen/xml`.
- `docs/conf.py` enables Breathe, MyST Markdown, MathJax, and the Read the Docs
  theme.
- `docs/index.rst` defines the published navigation.
- `.readthedocs.yaml` installs Doxygen and runs it before Sphinx.
- `docs/requirements.txt` pins the Python documentation dependencies.

When adding a public guide, place it under the appropriate `docs/` section and
add it to the toctree in `docs/index.rst`. Development notes should not be
presented as stable user instructions.

## Documenting Fortran

Use Doxygen comments immediately before a public routine:

```fortran
!> Compute a projected Hamiltonian.
!! @param[in]  h      Input Hamiltonian.
!! @param[out] hproj  Projected Hamiltonian.
subroutine project_hamiltonian(h, hproj)
```

Keep parameter descriptions consistent with the declaration and avoid
recording temporary debugging history in API comments.

## Troubleshooting

If API pages are empty, first check that Doxygen produced XML:

```sh
test -f docs/_build/doxygen/xml/index.xml
```

If a guide is absent from navigation, confirm that it is listed in the
toctree and that `myst-parser` is installed. If Read the Docs differs from a
local build, compare its Python and Doxygen versions with the configuration
above.
