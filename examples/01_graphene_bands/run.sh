#!/bin/sh
set -eu

example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
default_bin="$example_dir/../../lanczosKuboCode/bin/grabnes"
grabnes_bin=${GRABNES_BIN:-$default_bin}

if [ ! -x "$grabnes_bin" ]; then
    printf '%s\n' "GRABNES executable not found or not executable: $grabnes_bin" >&2
    printf '%s\n' "Compile GRABNES first or set GRABNES_BIN=/path/to/grabnes." >&2
    exit 1
fi

cd "$example_dir"
rm -f generate.bands job.out
"$grabnes_bin" Gendata.in > job.out

if [ ! -s generate.bands ]; then
    printf '%s\n' "GRABNES finished without creating generate.bands." >&2
    exit 1
fi

printf '%s\n' "Created $example_dir/generate.bands"
printf '%s\n' "Run 'python3 plot.py' to create graphene_bands.png."
