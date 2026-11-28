(development)=
# Development guide

The Fortran sources live in the `lanczosKuboCode/Src/` folder of the
[GRABNES code repository](https://github.com/gjung-group/grabnes). This
section collects the development documentation of that repository: the status
of every capability, the validation record, the conventions, and the design
notes behind recent changes.

## Code structure

The main interactions between the modules of the code are illustrated below.

![Main interactions in grabnes code](../assets/images/codeStructure.png)

The program and the modules making up the code are listed in the
[](../reference/index.md); the support libraries `math/` and `MIO/` sit
beside them in `lanczosKuboCode/Src/`.

## Testing and reproducibility

```sh
./tests/regression/smoke_test.sh      # build + example 01
./tests/regression/run_examples.sh    # build + examples 01 to 05 + further checks
python3 -m pytest tests/tapw         # TAPW linear-algebra conventions (Python only)
```

The same commands run automatically on GitHub for every push to the
development and main branches (`.github/workflows/regression.yml`).

The scripts build the solver in a temporary directory, never in the source
tree, run the example inputs, and compare the results numerically with the
reference data (bands within 2e-6 eV, DOS within 1e-9 relative). They also
compare the assembled Hamiltonian, entry by entry, with one constructed
independently in Python. Measured results for both compilers are in
[](../../docs/development/cluster-build-and-validation.md).

## What is supported

The support levels and the known limitations of the current code are on
[](../status.md); the full classification, capability by capability, is in
[](../../docs/development/functionality-status.md).

## Finding your way around the source

The quickest way to understand what an input key does is to look it up in
[](../../docs/user-guide/input-keys.md) — its "Read in" column names the source file
and routine — or to grep the sources:

```bash
grep -r "SuperCell ." lanczosKuboCode/Src/
```

## Building for development

Follow the [](../installation.md) instructions; the build is driven by
`make.sys` (see [](../../docs/development/building-documentation.md) for the documentation build).
The code is parallelized with OpenMP.

## Development documents

The documents below are copied unchanged from `docs/development/` of the code
repository.

Terms used throughout:

- **Implemented**: the code exists.
- **Tested**: it runs and passes defined checks in `tests/regression/`.
- **Independently validated**: results were checked against an analytic
  result or a second implementation.
- **Research / experimental**: present without sufficient validation for
  general use.

### Status and validation

What is tested and validated, and how the input keys were renamed.

| Document | Content |
| --- | --- |
| [](../../docs/development/functionality-status.md) | Every capability classified by the evidence that exists for it |
| [](../../docs/development/cluster-build-and-validation.md) | Build on Linux, test results for both compilers, defects fixed, independent Hamiltonian validation, parameters, performance |
| [](../../docs/development/input-key-renaming-proposal.md) | Table of the present and former names of the input keys; source of the renaming of October 2026 |

### Working on the code

Conventions and design notes for working on the solver sources.

| Document | Content |
| --- | --- |
| [](../../docs/development/coding-conventions.md) | Conventions of the solver sources (modules, routines, the `debug.log` trace, `Section.Name` input keys) and the rules added by the tests |
| [](../../docs/development/output-style.md) | Rules for everything the solver writes to `job.out`: module prefixes, section headers, debug output |
| [](../../docs/development/building-documentation.md) | How to build the Sphinx/Doxygen documentation |
| [](../../docs/development/soc-implementation.md), [](../../docs/development/tapw-chern-optimization.md), [](../../docs/development/known-issues/berry-curvature.md) | Notes on research functionality (spin-orbit terms, TAPW, Berry curvature) |
| [](../../docs/development/hamhopping-restructuring.md) | State and plan of the `HamHopping`/`HamOnSite` source restructuring |

### Project planning (internal notes)

The maintainers' notes for the first release; they are mirrored for completeness.

| Document | Content |
| --- | --- |
| [](../../docs/development/release-readiness.md) | Checklist for the first announced release: what is complete, what needs a maintainer decision, what is optional |
| [](../../docs/development/licensing-audit.md) | The license decision (GPL-3.0-or-later), the notices found, and the open points |
| [](../../docs/development/citation-checklist.md) | What `CITATION.cff` contains and what is deliberately absent |
| [](../../docs/development/software-paper-outline.md) | Outline of a possible software paper |
| [](../../docs/development/pre-announcement-changelog.md) | Chronological log of the consolidation work |
