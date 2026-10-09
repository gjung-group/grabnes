# TAPW mathematical tests

These tests exercise the linear-algebra conventions used by the TAPW code:
Fortran column ordering, label normalization, sparse two-step projection, and
orthonormal-subspace reconstruction. They do not execute the Fortran solver
and therefore complement, rather than replace, end-to-end reference examples.

From the repository root:

```sh
python3 -m pip install -r tests/requirements.txt
python3 -m pytest tests/tapw
```

All random data use explicit seeds. New tests should use numerical assertions,
state their tolerance, and avoid plots or interactive output.
