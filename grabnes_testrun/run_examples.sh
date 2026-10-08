#!/bin/sh
# Complete GRABNES example suite: build the canonical solver and reproduce the
# reference data of all four public examples. See README.md in this directory.
set -eu

. "$(dirname -- "$0")/lib/common.sh"
parse_args "$@"
harness_main

run_example 01_graphene_bands        generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 02_graphene_dos          generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"
run_example 03_twisted_bilayer_bands generate.bands    bands.dat bands "$bands_atol" "$bands_rtol"
run_example 04_twisted_bilayer_dos   generate.diag.DOS dos.dat   dos   "$dos_atol"   "$dos_rtol"

finish "the canonical solver was built and reproduces all four example references"
