# Output style of `job.out`

Everything the solver writes to standard output follows these rules. A new
message that breaks one of them is a bug, like a wrong number.

## 1. Every line names its module

Print with `call MIO_Print(message, 'module')`. The prefix is the Fortran
module the code lives in, in lower case:

`start`, `input`, `mio`, `cell`, `atoms`, `tbpar`, `neigh`, `ham`, `magf`,
`calc`, `kubo`, `diag`, `scf`, `gauss`, `hybrid`, `interface`, `parallel`,
`random`, `moire` (moireBLShift), and in the private build `semicl`, `hall`.

A routine name (`HamHopping:`) is never a prefix. `print *`,
`write(*,*)` and `write(6,...)` are not used outside `Src/MIO`. To print
numbers, write them into a character buffer with an explicit format
(`write(msg,'(a,f12.6,a)') 'E_F = ', ef, ' eV'`) or use `num2str`, then
print the buffer with `MIO_Print`.

Lines longer than 110 characters are wrapped by `MIO_Print`. Continuation
lines are indented by two spaces after the prefix. Keep a number and its unit
in one short message so they stay on one line.

## 2. Sections

Each stage of a run opens with a section header:

```
call MIO_Section('Neighbour search', 'neigh')
```

```

====== Neighbour search ============================================== [neigh]
```

(a blank line, then one line of `=` 78 characters wide with the title and the module).

The stages are: input and system description (`start`), cell and atoms,
tight-binding parameters, neighbour search, Hamiltonian, the calculation that
was selected (band structure, TAPW, Kubo/recursion DOS, Hall, spectral
function, semiclassics, ...), and the final summary (`mio`). Inside a
section, `call MIO_Subsection('TAPW basis', 'diag')` prints
`diag: --- TAPW basis ---`.

## 3. Wording and layout

- Plain sentences or `Name: value unit` items. Sentence case, no
  exclamation marks, no first person ("We use ...", "Let us ..."), no
  trailing period on an item.
- ASCII only: `Ang`, `deg`, `x`, `sqrt`, `hbar`, `moire`. No symbols,
  check marks or Greek letters, which several terminals and `grep` mangle.
- Units are always given. Energies printed in eV unless the line says g0.
- Numbers have a fixed format chosen for the quantity (not list-directed
  `*`), so the same line has the same layout in every run.
- Items that belong to the line above are indented by two spaces after the
  prefix:

  ```
  diag: TAPW basis
  diag:   N_G: 18
  diag:   Plane waves: 151
  ```

- A table has one header line and aligned rows.
- Input keys are written with their `Section.Name` form.

## 4. Warnings and errors

`MIO_Warning(message, 'module')` for anything the user should act on,
`MIO_Kill(message, 'module', 'Routine')` to stop. `MIO_Print` never prints
the words `ERROR` or `WARNING` (the regression tests search for them).

## 5. Debug output

Diagnostics that are useful only while developing or chasing a bug are
off by default:

```
call MIO_PrintDebug(message, 'diag')               ! module: DEBUG: message
call MIO_PrintDebug(message, 'diag', tapwDebug)    ! also when the module's own key is on
if (MIO_Debugging('diag') .or. tapwDebug) then     ! skip work that is only printed
```

They are switched on at run time, with no rebuild:

| key | effect |
|---|---|
| `Output.Debug .true.` | debug output of every module |
| `Output.DebugModules diag neigh` | debug output of the listed modules |
| a module's own key (`Diag.TAPWDebug`, `Diag.socDebug`) | that feature's debug output and files |

The compile-time `-DDEBUG` build keeps its routine entry/exit trace
(`MIO_Debug`) and the checks that are too slow for production.

## 6. Loops

Nothing is printed once per atom, bond or G-vector outside debug output.
Per k-point output is at most one progress line per k-point, or a progress
line every 10 % of a long loop.
