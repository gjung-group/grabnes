# Shared functions of the GRABNES test harness (POSIX sh; sourced, not executed).
#
# The harness never keeps a copy of the solver. It builds lanczosKuboCode/
# out of tree into a work directory, runs the inputs of examples/ there, and
# compares the results with examples/*/reference/.

harness_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$harness_dir/../.." && pwd)
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
        work_dir=$(mktemp -d "${TMPDIR:-/tmp}/grabnes_regression.XXXXXXXX") ||
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

# run_solver NAME INPUT OUTPUT
# Runs the solver on a copy of INPUT in $work_dir/run/NAME. Returns 1, after
# reporting and counting the failure, unless the run ends cleanly with OUTPUT.
run_solver() {
    name=$1; input=$2; output=$3
    run_dir="$work_dir/run/$name"

    printf '\n%s\n' "== $name"
    [ -f "$input" ] || die "missing input $input"
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$input" "$run_dir/Gendata.in"

    status=0
    # GRABNES_LAUNCHER is split on purpose (for example "mpirun -np 1").
    (cd "$run_dir" && ${GRABNES_LAUNCHER:-} "$grabnes_bin" Gendata.in > job.out 2> job.err) || status=$?
    if [ "$status" -ne 0 ]; then
        printf '%s\n' "  FAIL: the solver exited with status $status (see $run_dir)"
        grep 'ERROR' "$run_dir/job.out" | sed -n '1,3p' | sed 's/^/  /'
        sed -n '1,5p' "$run_dir/job.err" | sed 's/^/  /'
    elif ! grep -q '^0 errors, 0 warnings' "$run_dir/job.out"; then
        printf '%s\n' "  FAIL: the solver did not report '0 errors, 0 warnings' (see $run_dir/job.out)"
    elif [ ! -s "$run_dir/$output" ]; then
        printf '%s\n' "  FAIL: the solver did not write $output"
    else
        return 0
    fi
    failures=$((failures + 1))
    return 1
}

# run_example NAME OUTPUT REFERENCE KIND ATOL RTOL
run_example() {
    reference="$examples_dir/$1/reference/$3"
    kind=$4; atol=$5; rtol=$6
    [ -f "$reference" ] || die "missing reference $reference"
    run_solver "$1" "$examples_dir/$1/Gendata.in" "$2" || return 0
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

# Stochastic Kubo (Lanczos recursion) DOS of graphene against the exact DOS.
run_kubo_check() {
    case_dir="$harness_dir/kubo_graphene_dos"
    run_solver kubo_graphene_dos "$case_dir/Gendata.in" generate.DOS || return 0
    if python3 "$harness_dir/check_kubo_dos.py" --max-dev "$kubo_max_dev" --rms-dev "$kubo_rms_dev" \
            "$run_dir/generate.DOS" "$case_dir/reference/exact_dos.dat"; then
        printf '%s\n' "  PASS"
    else
        printf '%s\n' "  FAIL: the Kubo DOS deviates from the exact graphene DOS"
        failures=$((failures + 1))
    fi
}

# variant_input NAME SED_EXPRESSION...
# Writes $work_dir/inputs/NAME.in, a modified copy of the example-03 input.
variant_input() {
    variant_name=$1; shift
    mkdir -p "$work_dir/inputs"
    sed "$@" "$examples_dir/03_twisted_bilayer_bands/Gendata.in" > "$work_dir/inputs/$variant_name.in"
}

# As variant_input, starting from the pristine graphene input of example 01.
graphene_variant_input() {
    variant_name=$1; shift
    mkdir -p "$work_dir/inputs"
    sed "$@" "$examples_dir/01_graphene_bands/Gendata.in" > "$work_dir/inputs/$variant_name.in"
}

# run_hamiltonian_check NAME VERIFY_TABLES_OPTIONS...
# Runs $work_dir/inputs/NAME.in and compares the tables the solver writes with
# WriteDataFiles against the independent Hamiltonian of
# tools/hamiltonian/verify_tables.py (needs NumPy; skipped without it).
run_hamiltonian_check() {
    if ! python3 -c 'import numpy' >/dev/null 2>&1; then
        printf '\n%s\n%s\n' "== $1" "  SKIP: NumPy is not available"
        return 0
    fi
    check_name=$1; shift
    run_solver "$check_name" "$work_dir/inputs/$check_name.in" generate.s.mag || return 0
    verify_status=0
    python3 "$repo_root/tools/hamiltonian/verify_tables.py" "$run_dir" "$@" \
        > "$run_dir/verify_tables.log" 2>&1 || verify_status=$?
    grep -E 'neighbor entries:|missing from|without the reverse|max \|(H_solver - H_model|E_solver - E_model|H_ij)|largest (Hermiticity|eigenvalue)|E_model - E_solver|against the band file  |MISMATCH|Error|Traceback' \
        "$run_dir/verify_tables.log" | sed 's/^ */  /'
    if [ "$verify_status" -eq 0 ]; then
        printf '%s\n' "  PASS"
    else
        printf '%s\n' "  FAIL: see $run_dir/verify_tables.log"
        failures=$((failures + 1))
    fi
}

# Survey of every model switch (model_survey.py): each of its small cases must
# end as recorded (ok, refused, non-Hermitian, ...) and reproduce the stored
# fingerprint of its Hamiltonian or eigenvalues. Needs NumPy.
run_model_survey() {
    printf '\n%s\n' "== model_survey"
    if ! python3 -c 'import numpy' >/dev/null 2>&1; then
        printf '%s\n' "  SKIP: NumPy is not available"
        return 0
    fi
    survey_status=0
    python3 "$harness_dir/model_survey.py" --bin "$grabnes_bin" --work "$work_dir/model_survey" \
        --launcher "${GRABNES_LAUNCHER:-}" --check > "$work_dir/model_survey.log" 2>&1 || survey_status=$?
    grep -E '^(summary|DIFFERS|PASS|FAIL)' "$work_dir/model_survey.log" | sed 's/^/  /'
    if [ "$survey_status" -ne 0 ]; then
        printf '%s\n' "  FAIL: see $work_dir/model_survey.log"
        failures=$((failures + 1))
    fi
}

# Known-answer checks of features without a public example (check_physics.py):
# Landau levels of graphene, the band gap of hBN, sparse against dense
# diagonalization. Needs NumPy.
run_physics_checks() {
    if ! python3 -c 'import numpy' >/dev/null 2>&1; then
        printf '\n%s\n%s\n' "== physics checks" "  SKIP: NumPy is not available"
        return 0
    fi
    physics_status=0
    python3 "$harness_dir/check_physics.py" --bin "$grabnes_bin" --work "$work_dir/physics" \
        --launcher "${GRABNES_LAUNCHER:-}" > "$work_dir/physics.log" 2>&1 || physics_status=$?
    grep -vE '^(PASS|FAIL): ' "$work_dir/physics.log" | sed '/^$/N;/^\n$/D'
    failures=$((failures + $(grep -c '^  FAIL' "$work_dir/physics.log" || true)))
    if [ "$physics_status" -ne 0 ] && ! grep -q '^  FAIL' "$work_dir/physics.log"; then
        printf '%s\n' "  FAIL: see $work_dir/physics.log"
        failures=$((failures + 1))
    fi
}

# Symmetry check without stored reference: the four Dirac states of example 03
# at the moire K point (first k-point, bands 37-40) form two degenerate pairs.
# The historical neighbor search split them by 1e-4 eV.
run_degeneracy_check() {
    printf '\n%s\n' "== dirac_point_degeneracy"
    bands="$work_dir/run/03_twisted_bilayer_bands/generate.bands"
    if [ ! -s "$bands" ]; then
        printf '%s\n' "  FAIL: no bands from example 03"
        failures=$((failures + 1))
        return 0
    fi
    if python3 - "$bands" <<'PYEOF'
import sys
lines = open(sys.argv[1]).read().split("\n")
bands = int(lines[3].split()[0])
values = [float(x) for line in lines[4:] for x in line.split()][1:1 + bands]
e = values[bands // 2 - 2:bands // 2 + 2]
split = max(e[1] - e[0], e[3] - e[2])
print(f"  Dirac states at K (eV): {e[0]:.6f} {e[1]:.6f} {e[2]:.6f} {e[3]:.6f}   largest pair splitting: {split:.1e}")
sys.exit(0 if split <= 1e-6 else 1)
PYEOF
    then
        printf '%s\n' "  PASS"
    else
        printf '%s\n' "  FAIL: the Dirac states at K are not pairwise degenerate within 1e-6 eV"
        failures=$((failures + 1))
    fi
}

# The Hamiltonian of example 03, of its F2G2 variant, and of a cell smaller
# than the interlayer search radius, each against the independent model.
run_hamiltonian_checks() {
    variant_input hamiltonian_example03 -e 's/^WriteDataFiles .*/WriteDataFiles .true./'
    run_hamiltonian_check hamiltonian_example03 --intralayer="$nn_elements"
    variant_input hamiltonian_f2g2 -e 's/^WriteDataFiles .*/WriteDataFiles .true./' \
        -e 's/^TB.NeighLevels .*/TB.NeighLevels 5/'
    run_hamiltonian_check hamiltonian_f2g2 --intralayer="$f2g2_elements"
    variant_input hamiltonian_small_cell -e 's/^WriteDataFiles .*/WriteDataFiles .true./' \
        -e 's/^TB.NeighLevels .*/TB.NeighLevels 5/' \
        -e 's/^MoireCellParameters .*/MoireCellParameters 2 1 1 2/'
    run_hamiltonian_check hamiltonian_small_cell --intralayer="$f2g2_elements"
}

# The two-centre parameter vpppi0 depends on the intralayer model:
# - KoshinoIntralayer .true. (original Moon-Koshino model): default 2.7 eV, used
#   for the intralayer and the interlayer hopping;
# - F2G2-type intralayer models: default 3.5 eV, used in the interlayer pi term.
# Both defaults are checked against the independent model and in the log.
run_parameter_convention_checks() {
    variant_input koshino_intralayer_default -e 's/^WriteDataFiles .*/WriteDataFiles .true./' \
        -e 's/^TB.NeighLevels .*/TB.NeighLevels 5/' -e '$a KoshinoIntralayer .true.'
    run_hamiltonian_check koshino_intralayer_default --koshino-intralayer 5 --vpppi0 2.7
    expect_log_line koshino_intralayer_default 'Two-centre Vpppi0 = 2.7000 eV (intralayer and interlayer'
    expect_log_line hamiltonian_f2g2 'Two-centre Vpppi0 = 3.5000 eV (interlayer pi term only'
}

# expect_log_line RUN TEXT: the solver log of an earlier run must contain TEXT.
expect_log_line() {
    if [ -f "$work_dir/run/$1/job.out" ] && grep -qF "$2" "$work_dir/run/$1/job.out"; then
        printf '%s\n' "  PASS: $1 reports '$2'"
    else
        printf '%s\n' "  FAIL: $1 does not report '$2'"
        failures=$((failures + 1))
    fi
}

# Other systems. hBN monolayer: tables against the model t = 10.68 - 3.11 a with
# the default on-site energies. Twisted bulk (periodic along z) and the legacy
# NeighList routine: structure and consistency only.
run_other_system_checks() {
    graphene_variant_input hbn_monolayer -e 's/^WriteDataFiles .*/WriteDataFiles .true./' \
        -e 's/^TypeOfSystem .*/TypeOfSystem BoronNitride/'
    run_hamiltonian_check hbn_monolayer --g0 3.1 --intralayer=-3.0294 --onsite 3:3.09,4:-1.89 \
        --interlayer-cutoff 1.0 --periodic-z

    variant_input twisted_bulk -e 's/^WriteDataFiles .*/WriteDataFiles .true./' \
        -e 's/^CellHeight .*/CellHeight 6.68/' -e '$a Bulk .true.'
    run_hamiltonian_check twisted_bulk --structure-only 1 --periodic-z

    graphene_variant_input supercell_default -e 's/^SuperCell .*/SuperCell 4/' -e '/^nonBulkSmall/d'
    graphene_variant_input supercell_neighlist -e 's/^SuperCell .*/SuperCell 4/' -e '/^nonBulkSmall/d' \
        -e '$a Neigh.fastNNnotsquare .false.'
    run_solver supercell_default "$work_dir/inputs/supercell_default.in" generate.bands || return 0
    run_solver supercell_neighlist "$work_dir/inputs/supercell_neighlist.in" generate.bands || return 0
    if cmp -s "$work_dir/run/supercell_default/generate.bands" "$work_dir/run/supercell_neighlist/generate.bands"; then
        printf '%s\n' "  PASS: NeighList gives the bands of the default search for a 32-atom graphene cell"
    else
        printf '%s\n' "  FAIL: NeighList and the default search give different bands"
        failures=$((failures + 1))
    fi

    graphene_variant_input legacy_fastnn -e '$a Neigh.fastNN .true.'
    run_expect_error legacy_fastnn 'Neigh.fastNN and Neigh.fastNNnotsquareNotRectangle are not'
    graphene_variant_input legacy_notrectangle -e '$a Neigh.fastNNnotsquare .false.' \
        -e '$a Neigh.fastNNnotsquareNotRectangle .true.'
    run_expect_error legacy_notrectangle 'Neigh.fastNN and Neigh.fastNNnotsquareNotRectangle are not'
}

# run_expect_error NAME PATTERN: $work_dir/inputs/NAME.in must be refused with
# an error message matching PATTERN.
run_expect_error() {
    printf '\n%s\n' "== $1"
    run_dir="$work_dir/run/$1"
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$work_dir/inputs/$1.in" "$run_dir/Gendata.in"
    status=0
    (cd "$run_dir" && ${GRABNES_LAUNCHER:-} "$grabnes_bin" Gendata.in > job.out 2> job.err) || status=$?
    if [ "$status" -ne 0 ] && grep -q "$2" "$run_dir/job.out"; then
        printf '%s\n' "  PASS: refused with '$2'"
    else
        printf '%s\n' "  FAIL: expected the error '$2', got status $status (see $run_dir)"
        failures=$((failures + 1))
    fi
}

# TB.NeighLevels is the single shell control: values outside 1..8 are refused,
# and the deprecated Neigh.CutAtNN3 still means "three shells".
run_shell_control_checks() {
    variant_input neighlevels_0 -e 's/^TB.NeighLevels .*/TB.NeighLevels 0/'
    run_expect_error neighlevels_0 'TB.NeighLevels must be at least 1'
    variant_input neighlevels_9 -e 's/^TB.NeighLevels .*/TB.NeighLevels 9/'
    run_expect_error neighlevels_9 'TB.NeighLevels must be between 1 and 8'

    variant_input neighlevels_3 -e 's/^TB.NeighLevels .*/TB.NeighLevels 3/'
    variant_input legacy_cutatnn3 -e 's/^TB.NeighLevels .*/TB.NeighLevels 5/' -e '$a Neigh.CutAtNN3 .true.'
    run_solver neighlevels_3 "$work_dir/inputs/neighlevels_3.in" generate.bands || return 0
    run_solver legacy_cutatnn3 "$work_dir/inputs/legacy_cutatnn3.in" generate.bands || return 0
    if cmp -s "$work_dir/run/neighlevels_3/generate.bands" "$work_dir/run/legacy_cutatnn3/generate.bands" &&
       grep -q 'Neigh.CutAtNN3 is deprecated' "$run_dir/job.out"; then
        printf '%s\n' "  PASS: Neigh.CutAtNN3 with TB.NeighLevels 5 gives the bands of TB.NeighLevels 3"
    else
        printf '%s\n' "  FAIL: the legacy Neigh.CutAtNN3 input does not reproduce TB.NeighLevels 3"
        failures=$((failures + 1))
    fi
}

# A run with two MPI processes must stop with the documented error message
# instead of crashing: the MPI domain decomposition is disabled in the solver.
run_mpi_guard_check() {
    printf '\n%s\n' "== two_mpi_processes"
    mpirun_cmd=${GRABNES_MPIRUN:-mpirun}
    if ! command -v "${mpirun_cmd%% *}" >/dev/null 2>&1; then
        printf '%s\n' "  SKIP: '${mpirun_cmd%% *}' not found (set GRABNES_MPIRUN to the MPI launcher)"
        return 0
    fi
    run_dir="$work_dir/run/two_mpi_processes"
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$examples_dir/01_graphene_bands/Gendata.in" "$run_dir/"
    status=0
    (cd "$run_dir" && $mpirun_cmd -np 2 "$grabnes_bin" Gendata.in > job.out 2> job.err) || status=$?
    if [ "$status" -ne 0 ] && grep -q 'must be run with a single MPI' "$run_dir/job.out"; then
        printf '%s\n' "  PASS: the run was refused with the single-process error (status $status)"
    else
        printf '%s\n' "  FAIL: expected the single-process error, got status $status (see $run_dir)"
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

# The Kubo DOS uses one random-phase state on 180000 atoms. Six runs (three
# seeds, GNU and Intel random-number generators) gave rms deviations of
# 0.0027-0.0033 and maximum deviations of 0.008-0.014 from the exact DOS
# (peak 0.417). The limits are about twice that; a 2 % error in the hopping or
# a broadening in the wrong energy unit exceeds them.
kubo_max_dev=0.03
kubo_rms_dev=0.006

# Intralayer matrix elements (eV, one per neighbor shell) of the models checked
# by run_hamiltonian_checks; they restate the solver's default parameters and
# are passed to the independent Python model. verify_tables.py requires the
# matrix elements to agree within 1e-9 eV and the six-decimal band file within
# 1e-6 eV.
nn_elements=-2.9888
f2g2_elements=-2.9888,0.2354,-0.1877,0,0.0633
