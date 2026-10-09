# Development documentation

| Document | Content |
| --- | --- |
| [`release-readiness.md`](release-readiness.md) | Checklist for the first announced release: what is complete, what needs a maintainer decision, what is optional |
| [`functionality-status.md`](functionality-status.md) | Every capability classified by the evidence that exists for it |
| [`cluster-build-and-validation.md`](cluster-build-and-validation.md) | Build on Linux, test results for both compilers, defects fixed, independent Hamiltonian validation, parameters, performance |
| [`licensing-audit.md`](licensing-audit.md) | The license decision (GPL-3.0-or-later), the notices found, and the open points |
| [`citation-checklist.md`](citation-checklist.md) | What `CITATION.cff` contains and what is deliberately absent |
| [`software-paper-outline.md`](software-paper-outline.md) | Outline of a possible software paper |
| [`pre-announcement-changelog.md`](pre-announcement-changelog.md) | Chronological log of the consolidation work |
| [`building-documentation.md`](building-documentation.md) | How to build the Sphinx/Doxygen documentation |
| [`soc-implementation.md`](soc-implementation.md), [`tapw-chern-optimization.md`](tapw-chern-optimization.md), [`known-issues/`](known-issues/) | Notes on research functionality (spin-orbit terms, TAPW, Berry curvature) |

Terms used throughout:

- **Implemented**: the code exists.
- **Tested**: it runs and passes defined checks in `tests/regression/`.
- **Independently validated**: results were checked against an analytic
  result or a second implementation.
- **Research / experimental**: present without sufficient validation for
  general use.
