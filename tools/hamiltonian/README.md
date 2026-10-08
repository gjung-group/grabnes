# Hamiltonian verification

`verify_tables.py` checks the tight-binding Hamiltonian that GRABNES assembles
against an independent reconstruction in Python (NumPy required).

Run a band calculation with `WriteDataFiles .true.` in a scratch directory and
pass that directory to the script:

```sh
python3 tools/hamiltonian/verify_tables.py /path/to/run
```

It prints the intralayer shells with their matrix elements, compares the
interlayer elements with the two-center formula, reports how complete and how
symmetric the neighbor list is, and compares three reconstructions of the
bands with the solver's output. `--max-table-dev` and `--max-model-dev` turn
the last comparison into a pass/fail test; the regression suite uses them for
example 03 (`grabnes_testrun/run_examples.sh`).

The script assumes the band path of the public examples
(`K - Gamma - M - K'`) and the default two-center parameters; see `--help` for
the options. Results for the public examples are discussed in
[`docs/development/cluster-build-and-validation.md`](../../docs/development/cluster-build-and-validation.md).
