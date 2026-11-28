# Writing the input in YAML or JSON

The solver reads only the plain-text `Gendata.in` format described in
[](../../docs/user-guide/input-file.md). The code repository also ships a converter,
[`tools/input/gendata.py`](https://github.com/gjung-group/grabnes/blob/main/tools/input/gendata.py),
that turns a YAML or JSON description of a run into `Gendata.in` (`pack`) and
an existing `Gendata.in` back into YAML or JSON (`unpack`). It validates every
key against the list of keys the solver actually reads, so a misspelt key
becomes an error at conversion time instead of a key silently left at its
default.

```{note}
This page describes the converter on `main` of the code repository at commit
`2b5f86b`. The converter is a helper script: GRABNES itself still runs on
`Gendata.in`.
```

## Requirements

The converter is three Python scripts in `tools/input/` of a GRABNES checkout:

| File | Role |
| --- | --- |
| `gendata.py` | command-line interface: `pack`, `unpack`, `dump-schema` |
| `parser.py` | reads YAML/JSON data or `Gendata.in` text into one internal model |
| `gendata_schema.py` | generates the input schema from the Fortran sources |

It needs Python 3.10 or newer with `click` and `ruamel.yaml`; both are listed in
the code repository's `pyproject.toml`, so `uv sync` in the checkout installs
them. Run it from inside the checkout: by default the schema is generated from
`lanczosKuboCode/Src` of that checkout (see [](#the-schema)). No Fortran build is
needed.

```bash
python3 tools/input/gendata.py pack config.yaml -o Gendata.in
python3 tools/input/gendata.py unpack Gendata.in -o config.yaml
python3 tools/input/gendata.py dump-schema schema.json
```

## The YAML layout

Each `Section.Name` key becomes a key `name` inside a mapping `section`. Blocks
(`&begin ... &end`) become a list of rows. This is the input of
[example 01](../examples/getting-started.md) (graphene bands):

```yaml
run:
  description: Pristine graphene band structure
  typeOfSystem: Graphene
structure:
  cellSize: 1
  superCell: 1
tb:
  hopping: 3.1
  neighLevels: 1
diag:
  calc: true
calculate:
  bands: true
bands:
  numPoints: 20
  path:
    - [0.6666666666666666, 0.3333333333333333, 0.0]   # K
    - [0.0, 0.0, 0.0]                                 # Gamma
    - [0.5, 0.0, 0.0]                                 # M
    - [0.6666666666666666, 0.3333333333333333, 0.0]   # K'
```

`gendata.py pack` writes it as:

```
Run.Description Pristine graphene band structure
Run.TypeOfSystem Graphene

Structure.CellSize      1
Structure.SuperCell     1

TB.Hopping              3.1
TB.NeighLevels          1

Calculate.Bands         .true.

Diag.Calc               .true.

Bands.NumPoints         20
&begin Bands.Path 4
0.6666666666666666  0.3333333333333333  0.0   # K
0.0  0.0  0.0   # Gamma
0.5  0.0  0.0   # M
0.6666666666666666  0.3333333333333333  0.0   # K'
&end Bands.Path
```

Rules of the layout:

- **Names are matched as the solver matches them:** case is ignored, and `_` and
  `-` are dropped, so `typeOfSystem`, `TypeOfSystem` and `type_of_system` are the
  same key. `unpack` writes the section as in the source (`Run`, `TB`) and the
  key with a lower-case first letter (`typeOfSystem`).
- **Full names are accepted too.** A top-level entry whose value is not a
  mapping is taken as a full key name, so `Moire.TwistAngle: 1.05` and the
  section-less former names (`MoireTwistAngle: 1.05`) work.
- **Former names** (from before October 2026, see [](../../docs/user-guide/input-file.md)) are
  rewritten to the present name, with a warning.
- **Values keep their YAML type** and are checked against the type the solver
  expects: `true`/`false` (or `yes`, `no`, `t`, `f`, `.true.`, `.false.` as
  strings) for logical keys, numbers for integer and real keys, and lists of the
  right length for array keys such as `Diag.KGrid: [4, 4, 1]`. A string of
  arithmetic, as the solver accepts in `Gendata.in` (`2/3`, `0-0.01`), is
  evaluated; nothing else is executed.
- **Block rows** are lists of numbers; a block with a single row may be written
  as one flat list.
- **Comments** above a key, at the end of a key line, and at the end of a block
  row are carried into `Gendata.in`, and back when unpacking.

```{admonition} Observed limitation
:class: warning

A comment at the very top of a YAML file, before the first section, is not
carried into `Gendata.in` (tested at `2b5f86b`). Comments above later sections
and keys are kept.
```

## JSON

A JSON input is a flat object of full key names; blocks are lists of rows:

```json
{
 "Run.TypeOfSystem": "Graphene",
 "TB.Hopping": 3.1,
 "Diag.KGrid": [4, 4, 1]
}
```

`pack` decides the format from the extension (`.yaml`, `.yml` or `.json`); any
other extension is refused. `unpack --format json` writes this flat form. JSON
has no comments, so comments of the `.in` file are lost in that direction.

## What the conversion fixes for you

The quirks of the `.in` format listed in [](../../docs/user-guide/input-file.md) are handled when
writing it:

- `true`/`false` become `.true.`/`.false.`.
- **Negative numbers are written in the `0-` form**, for reals and integers
  alike: `-1.05` becomes `0-1.05` and `-5` becomes `0-5`. The solver reads
  integers with the same expression parser as reals, so a bare negative number
  stops a run with "Bad real number" in either case.
- Keys are grouped by section in a fixed order (`Run`, `Structure`, `TB`,
  `Moire`, ... , `Input`; unlisted sections come last), with the keys of a
  section sorted by name and its blocks after its keys.
- The row count of `&begin Name N` is computed from the rows.

For example, `moireTwistAngle: -1.05` and `kubo: {seedValue: -5}` give:

```
warning: moireTwistAngle: former name of Moire.TwistAngle, writing Moire.TwistAngle
Moire.TwistAngle        0-1.05

Kubo.SeedValue          0-5
```

## Errors and warnings

`pack` stops with exit status 1 and a one-line message, without writing a file,
when the input is wrong:

```
Error: moire.twistAngel: no key of this name in this version of the solver (did you mean Moire.TwistAngle?)
Error: structure.cellSize: expected an integer, got 'six'
Error: diag.kGrid: expected 3 numbers, got 2
```

Softer problems are printed as `warning:` lines on stderr and the file is still
written: a former name rewritten to the present one, or a key given both under
its present and a former name. For a repeated key the converter keeps the value
the solver would use: the present name wins over former names wherever it
stands, then the former names in the order of the solver's alias table, and
among identical spellings the first one.

`unpack` is lenient, since its input is a file that may already run: an unknown
key or block, a block whose `&begin Name N` count does not match its rows, and a
repeated key are warnings. It reads the file as the solver does, including
comments starting with `!` or `#` anywhere on a line, comment lines starting
with `!`, `#` or `*`, and lines continued with a trailing `\`. An `&end` without
`&begin`, a block that is never closed, or a key without a value is an error.

An unknown key that comes from `unpack` is refused by a later `pack` unless
`--allow-unknown` is given, which writes it unchanged with a warning. This
happens with real inputs: example 01's `Gendata.in` contains `WritePos`, which
the schema at `2b5f86b` does not list:

```
$ python3 tools/input/gendata.py unpack Gendata.in -o ex01.yaml
warning: WritePos: no key of this name in this version of the solver
wrote ex01.yaml: 16 keys and 1 blocks (1 warning(s))
$ python3 tools/input/gendata.py pack ex01.yaml
Error: writePos: no key of this name in this version of the solver
```

## Command reference

`pack CONFIG` converts a `.yaml`, `.yml` or `.json` file to `Gendata.in`.

| Option | Effect |
| --- | --- |
| `-o, --output FILE` | write to `FILE` (default: standard output) |
| `--schema FILE` | validate against a stored schema instead of the one generated from the source tree |
| `--allow-unknown` | write keys unknown to this version unchanged, with a warning |
| `--check` | only validate and print a summary (`g.yaml: 9 keys and 1 blocks, 0 warning(s)`); write nothing |

`unpack INFILE` converts a `Gendata.in` to YAML.

| Option | Effect |
| --- | --- |
| `-o, --output FILE` | write to `FILE` (default: standard output) |
| `--schema FILE` | use a stored schema |
| `--format yaml\|json` | output format (default `yaml`) |

`dump-schema [OUTPUT]` writes the schema of the source tree as JSON (to standard
output without `OUTPUT`). `python3 tools/input/gendata_schema.py [OUTPUT]` does
the same.

(the-schema)=
## The schema

Nothing in the converter is a hand-maintained list of keys. `gendata_schema.py`
extracts every key and block from the calls to `MIO_InputParameter` and
`MIO_InputBlock` in `lanczosKuboCode/Src`, the same extraction that generates
[](../../docs/user-guide/input-keys.md). Each version of the code therefore validates against exactly
the keys it reads. For every key the schema records:

- its present name and section;
- its type (`logical`, `integer`, `real`, `string` or `array`), inferred from
  the default expression in the source, and for arrays the item kind and length;
- its default expression(s) and the routines that read it;
- its former names.

It also lists the blocks (an indexed block such as `TB.Hoppings.` matches
`TB.Hoppings.2`, ...) and maps every former name to the present one. A stored
schema is stamped with the GRABNES version from `pyproject.toml`, the
`git describe` of the checkout and a hash of the Fortran sources:

```
$ python3 tools/input/gendata.py dump-schema schema.json
wrote schema.json: 753 keys, 6 blocks, grabnes 3.0b0 (2b5f86b)
```

A schema file can be archived next to an input, and passed with `--schema` to
check or convert an input against another version of the code. A schema
written in a different schema format version is refused.

## Tests

The code repository tests the converter in `tests/input/` (no Fortran build
needed):

```bash
uv run python -m pytest tests/input
```

Besides unit tests of the parsing and formatting rules above, the round-trip
tests convert every `examples/**/Gendata*.in` to `.in` and to YAML and back,
and check that the parsed values are unchanged.
