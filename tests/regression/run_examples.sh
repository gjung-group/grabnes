#!/bin/sh
# Complete GRABNES regression suite: build the canonical solver, reproduce the
# reference data of the public examples 01 to 05, and run the additional solver
# checks (independent Hamiltonian, neighbor-shell control, Kubo DOS, MPI
# process guard). See README.md in this directory.
set -eu

. "$(dirname -- "$0")/lib/common.sh"
parse_args "$@"
harness_main

run_example 01_graphene_bands        generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 02_graphene_dos          generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"
run_example 03_twisted_bilayer_bands generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 04_twisted_bilayer_dos   generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"

# Example 05 uses the plane-wave reduction, which must run with one OpenMP thread. Its levels come from a
# projected matrix, so builds with different linear-algebra libraries agree to a few 1e-6 eV only.
saved_omp_threads=${OMP_NUM_THREADS-}
OMP_NUM_THREADS=1; export OMP_NUM_THREADS
run_example 05_graphene_hbn_bands    generate.bands    bands.dat bands 5e-6 0
if [ -n "$saved_omp_threads" ]; then OMP_NUM_THREADS=$saved_omp_threads; else unset OMP_NUM_THREADS; fi

# Checks that are not public examples (inputs in this directory)
run_degeneracy_check
run_hamiltonian_checks
run_parameter_convention_checks
run_other_system_checks
run_shell_control_checks
run_kubo_check
run_mpi_guard_check
run_physics_checks
run_hermann_checks
run_model_survey
run_generated_files_check

finish "the canonical solver was built and passes all example and solver checks"
