#!/bin/sh
# Complete GRABNES regression suite: build the canonical solver, reproduce the
# reference data of all four public examples, and run the additional solver
# checks (Hamiltonian tables, Kubo DOS, MPI process guard). See README.md in this directory.
set -eu

. "$(dirname -- "$0")/lib/common.sh"
parse_args "$@"
harness_main

run_example 01_graphene_bands        generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 02_graphene_dos          generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"
run_example 03_twisted_bilayer_bands generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 04_twisted_bilayer_dos   generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"

# Checks that are not public examples (inputs in this directory)
run_hamiltonian_check
run_kubo_check
run_mpi_guard_check

finish "the canonical solver was built and passes all example and solver checks"
