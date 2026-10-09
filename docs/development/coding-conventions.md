# Coding conventions of the solver

These are the conventions of the original GRABNES sources (`kubo.F90`,
`kubosubs.F90`, `cell.F90`, `hybrid.F90`, `interface.F90`, `parallel.F90`,
`kuboarrays.f90`, and the `MIO` library), written down so that new code and the
clean-up of the later, larger files (`ham.F90`, `diag.F90`, `neigh.F90`,
`calc.F90`) follow one style. The second part lists rules that the tests of
October 2026 showed to be necessary.

## Conventions of the original code

**Files and modules**

- One module per file, named after the file. The module starts with `use mio`,
  `implicit none`, `PRIVATE`, and then one `public ::` line per exported name.
- Generic names are declared with an `interface` block
  (`KuboUpdate` for `KuboUpdate_d` and `KuboUpdate_z`).
- Module variables are few. Data shared between modules lives in the module
  that owns it (`atoms`, `cell`, `neigh`, `ham`) and is imported with
  `use module, only : names`.

**Routines**

- Names are `ModuleVerb` in mixed case: `KuboRecursion`, `CellGet`, `HybridGen`.
  Variables are lower camel case.
- A routine imports what only it needs with its own `use ..., only :` lines.
  The `only` keywords are aligned in one column.
- Declarations come in this order: dummy arguments, each with `intent`; a blank
  line; local variables. Reals are `real(dp)`, literals carry `_dp`, constants
  come from the `constants` module. Pointers are initialized with `=>NULL()`.
- A routine starts with the trace and timer hooks and ends with their
  counterparts:

  ```fortran
  #ifdef DEBUG
     call MIO_Debug('KuboRecursion',0)
  #endif /* DEBUG */
  #ifdef TIMER
     call MIO_TimerCount('kubo::Rec')
  #endif /* TIMER */
  ```

  Code that must always run is never placed inside the `DEBUG` block.
- Routines are short enough to read on a few screens (the largest original
  routine has 164 lines) and nest a few levels deep at most.

**The execution trace (`debug.log`)**

- A build with `-DDEBUG` writes the file `debug.log`. Every routine reports
  when it is entered and when it is left:

  ```fortran
  #ifdef DEBUG
     call MIO_Debug('KuboRecursion',0)     ! first executable statement
  #endif /* DEBUG */
     ...
  #ifdef DEBUG
     call MIO_Debug('KuboRecursion',1)     ! last executable statement
  #endif /* DEBUG */
  ```

  which gives the lines `DEBUG: Entering subroutine KuboRecursion` and
  `DEBUG: Finished subroutine KuboRecursion`. The file is flushed after every
  line, so after a crash its last line names the routine that was running. With
  MPI only the printing process writes.
- The name passed is the name of the routine, exactly. The exit call has the
  argument 1 and is placed before every `return` as well as at the end.
- The original files follow this throughout (`kubo.F90`, `gauss.F90`,
  `magf.F90`, `cell.F90`: every routine). The later files do not: 3 of 36
  routines in `ham.F90`, 16 of 119 in `diag.F90`, 2 of 17 in `neigh.F90`,
  none in `calc.F90`. A trace of those files is therefore incomplete, and the
  clean-up adds the two calls to every routine.
- The same pair exists for timing, `MIO_TimerCount('module::Routine')` and
  `MIO_TimerStop`, under `#ifdef TIMER`.

**Input keys**

- A key has the form `Section.Name`, both parts in mixed case with a capital
  first letter: `Diag.Calc`, `Kubo.Calc`, `Calculate.Bands`, `TB.Hopping`,
  `TB.NeighLevels`, `Neigh.LayerNeighbors`, `MagField.Integer`,
  `Bands.NumPoints`, `Spectral.WeiKu`. The section names the part of the code
  that reads the key (usually the module), so that related keys sort together
  and a key says where it belongs.
- 135 of the 722 keys follow this form, in 14 sections (`Diag` 47, `Spectral`
  20, `TB` 17, `Neigh` 12, `Calculate` 10, `MagField` 7, `Bands` 7, ...). The
  other 587 have no section, and 281 of them start with a lower-case letter
  (`realStrain`, `middleTwist`, `fourLayersSandwiched`, `useTAPW`).
- New keys follow the `Section.Name` form. The existing keys keep working: see
  "Renaming keys" below.

**Layout**

- Three spaces per indentation level; no tabs; no trailing blanks.
- Lines stay below 132 columns. A continued line ends with `&`, and the
  continuation is indented.
- `end subroutine Name`, `end do`, `end if` are written out.

**Services of the MIO library**

- Input: `call MIO_InputParameter('Key',variable,default)`; blocks with
  `MIO_InputFindBlock` and `MIO_InputBlock`.
- Output: `call MIO_Print('text','module')`. No bare `print *` or
  `write(*,*)` in finished code.
- Errors: `call MIO_Kill('what is wrong and what to do','module','Routine')`.
- Memory: `call MIO_Allocate(array,bounds,'name','module')` and
  `MIO_Deallocate`, so that the memory report is complete.
- Files: the `cl_file` type (`file%Open`, `file%GetUnit`, `file%Close`).

**OpenMP**

- Every parallel loop lists its `PRIVATE` and `REDUCTION` variables explicitly
  and is closed with `!$OMP END PARALLEL DO`.

**Units**

- Energies are in units of the nearest-neighbor hopping `g0` inside the solver;
  electron-volts appear only in the input and in the output.
- `hopp(j,i)` holds minus the hopping of neighbor entry `j` of atom `i`.

## Rules added by the tests of October 2026

Each rule answers a defect that was found (see
[`cluster-build-and-validation.md`](cluster-build-and-validation.md) and
`tests/regression/model_survey.py`).

1. **A model that cannot work stops.** If a combination of switches is not
   implemented, or needs data that are missing, the solver calls `MIO_Kill`
   with a message that names the switch and says what to do. It never continues
   with unset values.
2. **No parameter is used before it is set.** A parameter that only some
   branches set is marked as unset before them (see `t2KA` in `HamHopping`),
   and its use is checked.
3. **Every scalar assigned inside a parallel loop is private** unless all
   threads assign the same value. A shared scalar that takes different values
   in different iterations makes the result depend on the run.
4. **Loops over a per-thread range (`in1`, `in2`) must be valid on every
   thread**, including the threads that own no atoms.
5. **A key has one default.** The same key must not be read with different
   defaults in different routines; read it once and pass the value on.
   (`tools/input/list_input_keys.py` lists the present exceptions.)
6. **A switch that has no effect in the selected calculation is reported**
   with a warning, as is done for the spin terms outside TAPW.
7. **No commented-out code.** The history is in the version control system.
8. **Random numbers are seeded through `RandSeedFromInput`**, so that
   `setSeed .true.` reproduces a disordered Hamiltonian.
9. **Every new switch gets a case in `tests/regression/model_cases.py`**, and
   an intended change of a model is recorded with
   `model_survey.py --update-reference`.
10. **New behavior is opt-in.** A new switch defaults to the previous result;
    a default is changed only to replace a result that was wrong, and the
    change is stated in the log of the run.

## Renaming keys without breaking inputs

Bringing the 587 keys without a section to the `Section.Name` form must not
invalidate existing input files. The mechanism is in the input library
(`Src/MIO/input.F90`); the table of names (`Src/MIO/input_aliases.inc`) is
filled once the names in `input-key-renaming-proposal.md` are settled:

1. One table in the input library maps every new name to its former name
   (`Ham.RealStrain` to `realStrain`, `Stack.MiddleTwist` to `middleTwist`).
2. The solver asks for the new name. If the input file does not contain it,
   the library looks for the former name and uses its value. An input file with
   only old names therefore gives exactly the same run.
3. If both names are present, the new one is used and a warning is printed. If
   an old name was used, one line at the end of the run lists the old names
   found and their replacements; nothing stops.
4. The variables inside the solver are not renamed by this; only the strings
   passed to `MIO_InputParameter` change, which is a mechanical edit that
   `tools/input/list_input_keys.py` and the model survey check.

A pair is declared with `call InputAddAlias('Section.Name','formerName')` in
`input_aliases.inc`, or, for a trial, in the input file itself:

```
&begin Input.Aliases 1
Ham.Shells  TB.NeighLevels
&end Input.Aliases
```

A key is found under either name whichever of the two the solver asks for, so
the table can be filled before the strings in the sources are changed.

The former names stay accepted indefinitely unless a later release decides
otherwise and announces it.
