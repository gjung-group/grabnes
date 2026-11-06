#!/bin/sh
# Reproduce the reference data of the longer public examples (06 to 11): each takes minutes and
# example 07 needs 1.5 GB of memory, so they are not part of run_examples.sh.
#
#   GRABNES_BIN=/path/to/grabnes sh tests/regression/run_long_examples.sh [--work-dir DIR] [NAME ...]
#
# Every example is copied to the work directory and run there with its own run.sh; nothing is written
# into examples/. OMP_NUM_THREADS sets the threads of examples 06 and 07. Exit status 1 if a comparison fails.
#
# Examples 06 and 07 use the recursion method from one random vector with a fixed seed. The vector comes from the
# compiler's intrinsic random_number, so the same seed gives a different vector with another compiler, and the
# density of states then differs at the level of the noise of one vector (about 1e-3). Their references were made
# with GNU Fortran: with a binary from another compiler they are run and compared, but a difference is reported as
# not comparable instead of a failure. The compiler is read from the binary's Fortran runtime symbols
# (GRABNES_COMPILER=gnu|intel overrides it; when it cannot be told, the comparison stays strict).
set -eu

harness_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
examples_dir=$(CDPATH= cd -- "$harness_dir/../../examples" && pwd)
work_dir=""
names=""
while [ $# -gt 0 ]; do
    case $1 in
        --work-dir) work_dir=$2; shift ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) names="$names $1" ;;
    esac
    shift
done
[ -n "$names" ] || names="06_graphene_hbn_landau_levels 07_graphene_hbn_twist_series 08_magic_angle_twisted_bilayer 09_alternating_twist_trilayer 10_twisted_bilayer_spectral_function 11_alternating_twist_trilayer_spectral_function"
[ -n "$work_dir" ] || work_dir=$(mktemp -d)
mkdir -p "$work_dir"
work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
GRABNES_BIN=${GRABNES_BIN:-$examples_dir/../lanczosKuboCode/bin/grabnes}
case $GRABNES_BIN in /*) ;; *) GRABNES_BIN="$(pwd)/$GRABNES_BIN" ;; esac
export GRABNES_BIN

compiler=${GRABNES_COMPILER:-}
if [ -z "$compiler" ] && command -v strings > /dev/null 2>&1; then
    if strings "$GRABNES_BIN" | grep -q '^_gfortran_'; then compiler=gnu
    elif strings "$GRABNES_BIN" | grep -qE '^for_(write|read|open)'; then compiler=intel
    fi
fi
printf '%s\n' "Compiler of $GRABNES_BIN: ${compiler:-unknown}"

failures=0
not_comparable=0
for name in $names; do
    printf '\n%s\n' "== $name"
    run_dir="$work_dir/$name"
    rm -rf "$run_dir"
    cp -r "$examples_dir/$name" "$run_dir"
    if ! (cd "$run_dir" && ./run.sh > run.log 2>&1); then
        printf '%s\n' "  FAIL: run.sh ended with an error (see $run_dir/run.log)"
        failures=$((failures + 1)); continue
    fi
    case $name in
        06_*) (cd "$run_dir" && python3 plot.py --table result.dat > /dev/null)
              set -- --kind table --atol 1e-5 --rtol 0 "$run_dir/result.dat" "$run_dir/reference/dos_map.dat" ;;
        07_*) (cd "$run_dir" && python3 plot.py --table result.dat > /dev/null)
              set -- --kind table --atol 1e-5 --rtol 0 "$run_dir/result.dat" "$run_dir/reference/dos.dat" ;;
        10_*|11_*) (cd "$run_dir" && python3 plot.py --table result.dat > /dev/null)
              set -- --kind table --atol 2e-3 --rtol 0 "$run_dir/result.dat" "$run_dir/reference/summary.dat" ;;
        *)    set -- --kind bands --atol 5e-6 --rtol 0 "$run_dir/generate.bands" "$run_dir/reference/bands.dat" ;;
    esac
    rc=0
    python3 "$harness_dir/compare_output.py" "$@" || rc=$?
    if [ "$rc" -eq 0 ]; then
        printf '%s\n' "  PASS"
    else
        # only values outside the tolerance (status 1) can come from another random vector; a result that cannot
        # be read, has another shape or holds non-finite values fails with any compiler
        case $name in
            06_*|07_*) if [ "$rc" -eq 1 ] && [ -n "$compiler" ] && [ "$compiler" != gnu ]; then
                           printf '%s\n' "  NOT COMPARABLE: reference made with GNU Fortran; this $compiler build draws a" \
                               "  different random vector from the same seed (intrinsic random_number), see the header"
                           not_comparable=$((not_comparable + 1)); continue
                       fi ;;
        esac
        printf '%s\n' "  FAIL: the result does not reproduce the reference"
        failures=$((failures + 1))
    fi
done

printf '\n'
if [ "$failures" -eq 0 ] && [ "$not_comparable" -gt 0 ]; then
    printf '%s\n' "PASS: no failure; $not_comparable random-vector example(s) not comparable with this $compiler build, the others" \
        "reproduce their reference data (work directory: $work_dir)"
elif [ "$failures" -eq 0 ]; then
    printf '%s\n' "PASS: the long examples reproduce their reference data (work directory: $work_dir)"
else
    printf '%s\n' "FAIL: $failures long example(s) failed (work directory: $work_dir)"
    exit 1
fi
