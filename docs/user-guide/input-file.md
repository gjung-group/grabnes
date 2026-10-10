# The input file

The solver is started as `grabnes Gendata.in > job.out`. The input file is plain
text; this page describes its format and the messages the solver prints about
it. Every key is listed in [`input-keys.md`](input-keys.md).

## Format

```
# a comment
Run.TypeOfSystem      Graphene
Structure.CellSize    6
Structure.SuperCell   1
TB.NeighLevels        1
Diag.Calc             .true.

&begin Bands.Path 3
0.0        0.0        0.0
0.5        0.0        0.0
0.3333333  0.3333333  0.0
&end Bands.Path
```

- One key per line, followed by its value. Keys are not case sensitive, and
  `_` and `-` inside a key are ignored.
- Logical values: `.true.`, `true`, `t`, `yes`, `y` (and the corresponding
  forms of false).
- A block starts with `&begin Name N`, where `N` is the number of lines that
  follow, and ends with `&end Name`.
- **A negative number is written `0-0.01`.** A bare `-0.01` stops the run with
  "Bad real number". This also holds for the numbers inside a block.
- **If a key appears twice, the first value is used.** The solver prints a
  warning that names the key. When an input file is adapted from another one,
  edit the existing line instead of appending a new one.
- A key that is absent takes its default, which is listed in
  [`input-keys.md`](input-keys.md).

## Names of the keys

A key has the form `Section.Name` (`Structure.SuperCell`, `Strain.RealStrain`,
`Kubo.RecursionNumber`). The section says which part of the calculation the key
belongs to.

Until October 2026 most keys had no section (`SuperCell`, `realStrain`,
`RecursionNumber`). **Those former names keep working**: an input file may use
either name, and a file written before the change gives the same run. The
"Former name" column of [`input-keys.md`](input-keys.md) gives the
correspondence. If both names of a key are in the file, the present name is
used and the solver prints a warning.

## What the solver reports about the input

At the end of `job.out`:

- `input: note: N key(s) of the input file are given under a former name` —
  nothing to do. `Input.ListFormerNames .true.` lists them with their present
  names.
- `input: note: N key(s) of the input file were not read in this run`, followed
  by the keys. A key in this list is either **misspelt** or belongs to a part of
  the code that the selected calculation does not use (the Kubo keys in a band
  calculation, for instance). Look through the list after the first run of a new
  input: a misspelt key is otherwise taken at its default without any other
  sign.

During the run:

- `input: WARNING: the key "..." appears N times` — see above.
- `ERROR` followed by `*** PROGRAM ABORTED` — the run stopped; the message says
  which switch or value is the reason. Combinations of model switches that are
  not implemented, or that would use values that were never set, stop in this
  way instead of producing numbers.
- `WARNING` lines do not stop the run. The ones to take seriously are those
  about a hopping table that is not Hermitian and about neighbor pairs close to
  the search radius.

## Things that change a run more than one expects

- `MagField.InitialInteger` and `MagField.FinalInteger` start a sweep over
  magnetic fields, one calculation per value, each in its own directory. For a
  single field use `MagField.Integer`; `MagField.Integer 0` is zero field.
- `Structure.SuperCell N` repeats the cell `N x N`; time and memory grow with
  the number of atoms.
- `TB.Hopping` sets the energy unit of the whole Hamiltonian and differs between
  families of models; do not carry it from one kind of input to another.
- With the plane-wave reduction (`TAPW.Use .true.`), `Diag.N_G` must be set:
  its default gives a very large basis.

## Execution trace

A solver built with `-DDEBUG` (add it to `FCFLAGS` in `make.sys`) writes
`debug.log`, one line when a routine is entered and one when it is left, and
prints additional diagnostic output. After a crash the last line of `debug.log`
names the routine that was running.
