# Citation metadata checklist

A `CITATION.cff` file was **not** created, because the information it must
contain cannot be established from the repository. Nothing below is to be read
as a statement of authorship.

## What the repository shows

| Item | Evidence in the repository |
| --- | --- |
| Name | GRABNES (GRAphene and Boron Nitride Electronic Structure) |
| Repository | `https://github.com/gjung-group/grabnes` |
| Version string | `3.0b` in `lanczosKuboCode/version.info`; `3.0b0` in `pyproject.toml` |
| Named in source headers | Rafael Martinez-Gordillo (copyright 2012, MIO and math libraries; "original" author in `ham.F90`); "Nicolas (modifications)" and "Jinwoo (current version)" in `ham.F90`, first names only |
| Sphinx metadata | author "Jeil Jung Group"; copyright "2022, Nicolas Leconte" |
| Git history | commits by the accounts "George Jung" and "Deepanshu Aggarwal", and the present development branch |
| Publications named in the repository | Jung, Raoux, Qiao, and MacDonald, Phys. Rev. B 89, 205414 (2014), for the effective graphene/hBN model (`ham.F90`, tutorial notebook); A. Cresti, Phys. Rev. B 103, 045402 (2021), for the Peierls phase (tutorial notebook); Moon and Koshino, Phys. Rev. B 85, 195458 (2012), for the interlayer two-center form (example README) |

## What has to be supplied by the maintainers

- [ ] The complete author list, with full names, order, affiliations, and
      ORCID identifiers where available.
- [ ] Whether earlier contributors named in the sources are authors of the
      release or are acknowledged.
- [ ] The license identifier (see `licensing-audit.md`); `CITATION.cff`
      requires one.
- [ ] The version number and date of the first release.
- [ ] A DOI for the release (for example from an archive such as Zenodo),
      once a release exists.
- [ ] The reference to be cited for the code itself (a software paper or a
      methods paper), once it exists.
- [ ] The list of publications that should be cited for individual models
      and methods, verified against the publications.

## Template

To be completed and placed at the top level as `CITATION.cff` when the items
above are known. Every value in angle brackets is a placeholder.

```yaml
cff-version: 1.2.0
message: "If you use GRABNES, please cite it as below."
title: "GRABNES: GRAphene and Boron Nitride Electronic Structure"
type: software
authors:
  - family-names: "<family name>"
    given-names: "<given names>"
    affiliation: "<affiliation>"
    orcid: "<https://orcid.org/...>"
repository-code: "https://github.com/gjung-group/grabnes"
license: "<SPDX identifier>"
version: "<version>"
date-released: "<YYYY-MM-DD>"
doi: "<DOI>"
```
