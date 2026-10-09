#!/bin/sh
# Minimal GRABNES check: build the canonical solver and reproduce the pristine
# graphene band structure. See README.md in this directory.
set -eu

. "$(dirname -- "$0")/lib/common.sh"
parse_args "$@"
harness_main

run_example 01_graphene_bands generate.bands bands.dat bands "$bands_atol" "$bands_rtol"

finish "the canonical solver was built and reproduces the graphene reference bands"
