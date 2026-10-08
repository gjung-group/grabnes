# Shared functions of the GRABNES test harness (POSIX sh; sourced, not executed).
#
# The harness never keeps a copy of the solver. It builds lanczosKuboCode/
# out of tree into a work directory, runs the inputs of examples/ there, and
# compares the results with examples/*/reference/.

harness_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$harness_dir/.." && pwd)
solver_dir="$repo_root/lanczosKuboCode"
examples_dir="$repo_root/examples"

keep_work=0
make_sys=${GRABNES_MAKE_SYS:-}
exact=0
failures=0

die() {
    printf '%s\n' "ERROR: $*" >&2
    # Do not leave an unused temporary work directory behind.
    if [ "${keep_work:-1}" -eq 0 ] && [ -n "${work_dir:-}" ]; then
        rmdir "$work_dir" 2>/dev/null || true
    fi
    exit 2
}

usage() {
    cat <<USAGE
Usage: $(basename "$0") [options]

Builds the canonical solver (lanczosKuboCode/) out of tree and checks example
results against their reference data.

Options:
  --debug            build with config/gfortran.debug.make.sys (runtime checks)
  --make-sys FILE    build configuration to use
  --work-dir DIR     build and run in DIR and keep it (allows incremental
                     rebuilds); the default is a new temporary directory that
                     is removed when every check passes
  --keep             keep the temporary work directory
  --exact            additionally require byte-identical band files
  -h, --help         show this help

Environment:
  GRABNES_MAKE_SYS   same as --make-sys
  GRABNES_TEST_DIR   same as --work-dir
  GRABNES_BIN        test this existing executable instead of building one
  GRABNES_LAUNCHER   command prefix used to start the executable, for example
                     "mpirun -np 1" or "srun -n 1" (default: none)

Build configuration, in order of precedence: --debug or --make-sys,
GRABNES_MAKE_SYS, lanczosKuboCode/make.sys if it exists, and otherwise
lanczosKuboCode/make.sys.example.
USAGE
}

parse_args() {
    work_dir=${GRABNES_TEST_DIR:-}
    while [ $# -gt 0 ]; do
        case $1 in
            --debug) make_sys="$harness_dir/config/gfortran.debug.make.sys" ;;
            --make-sys) [ $# -ge 2 ] || die "--make-sys needs a file"; make_sys=$2; shift ;;
            --work-dir) [ $# -ge 2 ] || die "--work-dir needs a directory"; work_dir=$2; shift ;;
            --keep) keep_work=1 ;;
            --exact) exact=1 ;;
            -h|--help) usage; exit 0 ;;
            *) usage >&2; die "unknown option: $1" ;;
        esac
        shift
    done
}

need_command() {
    command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH${2:+ ($2)}"
}

# Value of a variable defined in the selected make.sys.
make_sys_value() {
    printf 'include %s\nharness-print:\n\t@echo $(%s)\n' "$make_sys" "$1" |
        make -s -f - harness-print 2>/dev/null
}

prepare_work_dir() {
    if [ -n "$work_dir" ]; then
        mkdir -p "$work_dir" || die "cannot create work directory $work_dir"
        work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
        keep_work=1
    else
        work_dir=$(mktemp -d "${TMPDIR:-/tmp}/grabnes_testrun.XXXXXXXX") ||
            die "cannot create a temporary work directory (set GRABNES_TEST_DIR)"
    fi
    case $work_dir/ in
        "$solver_dir"/*|"$examples_dir"/*)
            die "the work directory must not be inside lanczosKuboCode/ or examples/" ;;
    esac
    printf '%s\n' "Work directory: $work_dir"
}

# Report the state of the checkout so that a result can be traced to a commit.
describe_source() {
    if command -v git >/dev/null 2>&1 && git -C "$repo_root" rev-parse HEAD >/dev/null 2>&1; then
        commit=$(git -C "$repo_root" rev-parse HEAD)
        if [ -n "$(git -C "$repo_root" status --porcelain -- lanczosKuboCode examples)" ]; then
            commit="$commit (with uncommitted changes in lanczosKuboCode/ or examples/)"
        fi
        printf '%s\n' "Source commit:  $commit"
    else
        printf '%s\n' "Source commit:  unknown (not a Git checkout)"
    fi
}

build_solver() {
    need_command make "GNU Make is required"
    make --version 2>/dev/null | grep -q 'GNU Make' || die "GNU Make is required"

    if [ -z "$make_sys" ]; then
        if [ -f "$solver_dir/make.sys" ]; then
            make_sys="$solver_dir/make.sys"
        else
            make_sys="$solver_dir/make.sys.example"
        fi
    fi
    [ -f "$make_sys" ] || die "build configuration not found: $make_sys"
    make_sys=$(CDPATH= cd -- "$(dirname -- "$make_sys")" && pwd)/$(basename -- "$make_sys")

    f90=$(make_sys_value F90)
    f90_serial=$(make_sys_value F90_SERIAL)
    [ -n "$f90" ] || die "F90 is not defined in $make_sys"
    need_command "${f90%% *}" "compiler named by F90 in $make_sys"
    [ -z "$f90_serial" ] || need_command "${f90_serial%% *}" "compiler named by F90_SERIAL in $make_sys"

    printf '%s\n' "Configuration:  $make_sys"
    printf '%s\n' "Compiler:       $f90 -> $("${f90%% *}" --version 2>/dev/null | sed -n 1p)"
    printf '%s\n' "FCFLAGS:        $(make_sys_value FCFLAGS)"
    printf '%s\n' "LIBS, LDFLAGS:  $(make_sys_value LIBS) $(make_sys_value LDFLAGS)"

    printf '%s\n' "Building lanczosKuboCode/ (log: $work_dir/build.log) ..."
    if ! make -C "$solver_dir" MAKE_SYS="$make_sys" \
            BUILD_DIR="$work_dir/build" BIN_DIR="$work_dir/bin" \
            > "$work_dir/build.log" 2>&1; then
        printf '%s\n' "---- last lines of build.log ----" >&2
        tail -n 30 "$work_dir/build.log" >&2
        printf '%s\n' "FAIL: compilation failed. The usual causes are a missing MPI wrapper," >&2
        printf '%s\n' "      BLAS/LAPACK or ARPACK library; adapt a copy of make.sys.example." >&2
        exit 1
    fi
    grabnes_bin="$work_dir/bin/grabnes"
    [ -x "$grabnes_bin" ] || die "the build finished without creating $grabnes_bin"
}

select_executable() {
    describe_source
    if [ -n "${GRABNES_BIN:-}" ]; then
        grabnes_bin=$GRABNES_BIN
        [ -x "$grabnes_bin" ] || die "GRABNES_BIN is not an executable: $grabnes_bin"
        printf '%s\n' "Testing the existing executable (no build): $grabnes_bin"
    else
        build_solver
    fi
    if command -v md5sum >/dev/null 2>&1; then
        printf '%s\n' "Executable md5: $(md5sum "$grabnes_bin" | cut -d' ' -f1)"
    fi
}

# run_example NAME OUTPUT REFERENCE KIND ATOL RTOL
run_example() {
    name=$1; output=$2; reference="$examples_dir/$1/reference/$3"
    kind=$4; atol=$5; rtol=$6
    run_dir="$work_dir/run/$name"

    printf '\n%s\n' "== $name"
    [ -f "$examples_dir/$name/Gendata.in" ] || die "missing input $examples_dir/$name/Gendata.in"
    [ -f "$reference" ] || die "missing reference $reference"
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$examples_dir/$name/Gendata.in" "$run_dir/"

    status=0
    # GRABNES_LAUNCHER is split on purpose (for example "mpirun -np 1").
    (cd "$run_dir" && ${GRABNES_LAUNCHER:-} "$grabnes_bin" Gendata.in > job.out 2> job.err) || status=$?
    if [ "$status" -ne 0 ]; then
        printf '%s\n' "  FAIL: the solver exited with status $status (see $run_dir)"
        grep 'ERROR' "$run_dir/job.out" | sed -n '1,3p;s/^/  /'
        sed -n '1,5p' "$run_dir/job.err" | sed 's/^/  /'
        failures=$((failures + 1))
        return 0
    fi
    if ! grep -q '^0 errors, 0 warnings' "$run_dir/job.out"; then
        printf '%s\n' "  FAIL: the solver did not report '0 errors, 0 warnings' (see $run_dir/job.out)"
        failures=$((failures + 1))
        return 0
    fi
    if [ ! -s "$run_dir/$output" ]; then
        printf '%s\n' "  FAIL: the solver did not write $output"
        failures=$((failures + 1))
        return 0
    fi
    set -- --kind "$kind" --atol "$atol" --rtol "$rtol"
    if [ "$exact" -eq 1 ] && [ "$kind" = bands ]; then
        set -- "$@" --exact
    fi
    if python3 "$harness_dir/compare_output.py" "$@" "$run_dir/$output" "$reference"; then
        printf '%s\n' "  PASS"
    else
        printf '%s\n' "  FAIL: $output does not reproduce $reference"
        failures=$((failures + 1))
    fi
}

finish() {
    printf '\n'
    if [ "$failures" -ne 0 ]; then
        printf '%s\n' "FAIL: $failures check(s) failed. Logs and outputs are kept in $work_dir"
        exit 1
    fi
    if [ "$keep_work" -eq 1 ]; then
        printf '%s\n' "PASS: $1 (work directory kept: $work_dir)"
    else
        rm -rf "$work_dir"
        printf '%s\n' "PASS: $1"
    fi
}

harness_main() {
    need_command python3 "needed to compare results with the reference data"
    prepare_work_dir
    select_executable
}

# Tolerances.
#
# Band files are written with six decimals; they have been byte-identical for
# every compiler tested so far, and 2e-6 only allows a last-digit rounding flip.
# The DOS files carry 17 significant digits, so agreement is limited by the
# BLAS/LAPACK and math-library rounding (observed: 1e-14 relative).
bands_atol=2e-6
bands_rtol=0
dos_atol=1e-12
dos_rtol=1e-9
