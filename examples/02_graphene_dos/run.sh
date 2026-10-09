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
rm -f generate.diag.DOS job.out file.13 file.14
"$grabnes_bin" Gendata.in > job.out

if [ ! -s generate.diag.DOS ]; then
    printf '%s\n' "GRABNES finished without creating generate.diag.DOS." >&2
    exit 1
fi

printf '%s\n' "Created $example_dir/generate.diag.DOS"
printf '%s\n' "Run 'python3 plot.py' to create graphene_dos.png."
