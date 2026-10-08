# Citation metadata

`CITATION.cff` at the top level is the citation record of GRABNES. It follows
the Citation File Format 1.2.0 and was validated against the official schema.

## What it contains

| Field | Value | Source |
| --- | --- | --- |
| Title | GRABNES: Graphene and Boron Nitride Electronic Structure | maintainer |
| Authors, in order | Jeil Jung; Rafael Martinez-Gordillo; Nicolas Leconte | maintainer |
| Contact | Nicolas Leconte (current developer and maintainer) | maintainer |
| Repository | `https://github.com/gjung-group/grabnes` | repository |
| License | `GPL-3.0-or-later` | `licensing-audit.md` |

Authorship records who wrote the software and in which order they are to be
cited. Rafael Martinez-Gordillo wrote the original program and its libraries;
the code was developed further in the group of Jeil Jung; Nicolas Leconte is
the sole current developer and the maintainer. The file makes no statement
about copyright ownership.

## Deliberately absent

These fields were left out because no reliable value exists yet. They must
not be filled with guesses.

- [ ] ORCID identifiers, affiliations, and e-mail addresses of the authors.
- [ ] `version` and `date-released`: no release has been made. (The sources
      carry the internal version string 3.0b.)
- [ ] `doi`: requires an archived release.
- [ ] `preferred-citation`: to be added when a software paper is published.
      None exists; a manuscript is only outlined in
      `software-paper-outline.md`.
- [ ] References to the publications behind individual models, verified
      against the publications.

## How to cite until then

Authors, title, repository URL, and the exact Git commit used, as shown in the
README.

## Validation

```sh
python3 - <<'PY'
import json, urllib.request, yaml, jsonschema
schema = json.load(urllib.request.urlopen(
    "https://raw.githubusercontent.com/citation-file-format/citation-file-format/1.2.0/schema.json"))
jsonschema.validate(yaml.safe_load(open("CITATION.cff")), schema)
print("valid")
PY
```
