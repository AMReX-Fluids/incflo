# Regression testing

incflo uses the [AMReX regression testing tools](https://github.com/AMReX-Codes/regression_testing)
to build problems, run them, compare plotfiles with stored benchmark solutions,
and generate HTML reports.

The upstream nightly suite spans all four build directories:

| Directory | Dimension | Embedded boundaries |
| --- | --- | --- |
| `test_2d` | 2 | Enabled |
| `test_3d` | 3 | Enabled |
| `test_no_eb_2d` | 2 | Disabled |
| `test_no_eb_3d` | 3 | Disabled |

Problem inputs and GNUmakefiles stay in these directories. Each configuration
section selects one through `buildDir`, relative to the incflo source checkout,
and selects an input within it through `inputFile`.

## Example and upstream suite

[incflo-tests.ini.example](incflo-tests.ini.example) is a starter subset with one
case per build configuration. Its case settings come from the
[November 5, 2026 nightly configuration](https://ccse.lbl.gov/pub/RegressionTesting1/incflo/2026-11-05/incflo-tests.ini).
The example uses distinct case names, placeholder paths, and disabled email.
Its benchmark names therefore differ from those of the upstream cases.

The configuration used by the upstream runner defines nightly coverage; this
example does not reproduce the full suite. Consult the
[upstream regression reports](https://ccse.lbl.gov/pub/RegressionTesting1/incflo/)
and the configuration archived with the relevant run when checking coverage or
reproducing a failure. An archived configuration also needs local paths and
reporting settings adjusted before use.

## Set up a scratch area

Use dedicated scratch clones: the runner can update repositories and check out
configured branches or commits. You need Git, Python with the dependencies
listed by the regression testing tools, GNU Make, a suitable C++ compiler,
and MPI for the example configuration.

```sh
export REGTEST_SCRATCH=/absolute/path/to/regtest-scratch
mkdir -p "${REGTEST_SCRATCH}/test_data/incflo"
git clone https://github.com/AMReX-Codes/amrex.git "${REGTEST_SCRATCH}/amrex"
git clone https://github.com/AMReX-Fluids/AMReX-Hydro.git "${REGTEST_SCRATCH}/AMReX-Hydro"
git clone https://github.com/AMReX-Codes/regression_testing.git "${REGTEST_SCRATCH}/regression_testing"
git clone https://github.com/AMReX-Fluids/incflo.git "${REGTEST_SCRATCH}/incflo"
cp "${REGTEST_SCRATCH}/incflo/regression/incflo-tests.ini.example" \
   "${REGTEST_SCRATCH}/incflo-tests.ini"
```

Edit the copied `incflo-tests.ini` outside the checkout:

- Replace every `/absolute/path/to/regtest-scratch` with the actual absolute
  scratch path. Shell variables are not expanded in these INI values.
- Set `dir` and `branch` (or `hash`) for AMReX, incflo, and AMReX-Hydro.
- Set `testTopDir` and `webTopDir` to your output and report locations.
- Adjust `COMP`, `numMakeJobs`, `MPIcommand`, and each case's
  `numprocs` to suit your machine.

The sibling clone layout matches the GNUmakefile defaults. For a different
layout, set `AMREX_HOME` and `AMREX_HYDRO_HOME` to the dependency locations.
The `[extra-AMReXHydro]` section lets the runner manage that dependency too.
Email and the report's parent-page link are disabled in the example; configure
these explicitly if needed for a hosted runner.

## Create benchmarks and run comparisons

Create initial reference solutions before running comparisons:

```sh
cd "${REGTEST_SCRATCH}"
python regression_testing/regtest.py --make_benchmarks "Initial local references" incflo-tests.ini
```

Benchmark creation establishes the expected answers; review the results before
accepting them. Repeating this command replaces existing benchmarks for the
selected cases. To compare subsequent runs against those references:

```sh
python regression_testing/regtest.py incflo-tests.ini
```

With the example paths, reports appear at
`${REGTEST_SCRATCH}/test_data/incflo/web/index.html`.
By default the runner updates the configured repositories. To test the current
state of the scratch clones, check out the desired revisions there and run:

```sh
python regression_testing/regtest.py --no_update All incflo-tests.ini
```

This tests the configured scratch clones. Copy or commit your intended changes
into those clones before running. Keep baseline and comparison dependency
revisions consistent when investigating a change. Consult the tool's help for
additional options:

```sh
python regression_testing/regtest.py -h
```

## Add or change a regression case

Add the input to the appropriate `test*` directory. Define a uniquely named
section in your runner configuration with matching `buildDir`, `inputFile`, and
`dim`; set process counts, runtime overrides, comparison tolerances, and restart
settings as needed. Build and run it, inspect its output, then create a reviewed
benchmark.

For upstream nightly inclusion, coordinate the runner configuration and
benchmark update with the incflo maintainers through the contribution review.
Adding an input or editing this example alone does not register a nightly test.
The upstream runner's maintainers manage its deployed configuration and stored
benchmarks. Update this guide or the starter subset when their documented setup
needs to change; keep machine-specific paths, recipients, and benchmark data
outside the repository.
