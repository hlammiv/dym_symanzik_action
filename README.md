# dym_symanzik_action

Discrete Yang–Mills lattice gauge theory with a **Symanzik-improved action**.

Monte Carlo simulation of pure gauge theory in which the gauge group is a *finite
subgroup* of SU(N) (primarily **S(1080) ⊂ SU(3)**). The action combines the
standard Wilson plaquette term (`beta1`) with a higher-order rectangular-loop
(Symanzik improvement) term (`beta2`):

    S = beta0 + beta1 * sum_plaquettes ReTr(U) + beta2 * sum_rectangles ReTr(U)

Observables measured each configuration: action/plaquette density, Polyakov
loops, and an `Nx × Nt` grid of Wilson loops. The goal is to map the
`(beta1, beta2)` phase diagram and study how well the discrete group approximates
continuous SU(3).

## Layout

    *.c, *.cpp, *.h     Simulation core (C/C++)
    Makefile            Build (Makefile_mac for macOS)
    groups/             Group multiplication tables + generators (see groups/README.md)
    verify_group.c      Standalone format checker for group files
    scripts/            Python analysis/plotting
    results/            Figures (PNG/PDF) and graph outputs      [git-ignored]
    data/               Raw Monte Carlo run logs                 [git-ignored]
    evermore.sh         Run driver
    convert.sh          Wraps kentucky2nersc
    kentucky2nersc      Precompiled config converter (SU(3) only; no source in tree)

### Core sources

    dym-mod-metro.cpp            Main MC driver (Metropolis + OpenMP; measures observables)
    dym-mod-metro-savecfg.cpp    Variant that dumps configurations in matrix form (SU(3) only)
    dym-mod-metroOG.cpp          Earlier reference version
    dym-mod-mod-metro.cpp        Variant
    rect_test.cpp /
      rect_test_computation.cpp  Rectangle-loop (improvement term) tests
    group.c / group.h            Group-table loading: mult table, Re/Im traces, inverses
    lattice.c / lattice.h        Lattice geometry; plaquette, Polyakov, Wilson loops
    wilson_flow.c                Wilson-flow scale setting (naive SU(3))
    timer.cpp / timer.h          Per-step timing

## Build

    make dym-mod-metro

Builds with `-O3 -fopenmp`. Requires `g++` and OpenMP.

The largest supported group order is `PMAX` in `group.h` (currently **5040**, to
fit the SU(4) subgroups). The `mult` table is `PMAX²` ints but demand-paged, so a
large `PMAX` costs nothing for small-group runs; bump it only to admit a bigger
group. See `groups/README.md` for the group-file format and the full catalogue.

## Run

    ./evermore.sh {beta1} {beta2} {D} {Nt} {Nx} {group_file} {output_folder/}

`group_file` is a filename under `groups/` (e.g. `mys1080-v4`). `output_folder`
must end with `/`. Example:

    ./evermore.sh 1.5 0.1 4 2 2 mys1080-v4 data/

Measurement lines in the output log are tagged `GMES:`.

The binary takes optional trailing args to override the sweep counts without
recompiling (defaults `1000 / 1000 / 0`):

    ./dym-mod-metro <group> D Nt Nx beta0 beta1 beta2 seed [K] [N] [Ntherm]

where `K` = decorrelation sweeps between measurements, `N` = number of
measurements, `Ntherm` = thermalization sweeps. `dym-mod-metro-savecfg` (which
writes SU(3) configs in Kentucky format for Wilson-flow scale setting, and needs
a group file with the defining-rep matrices appended, e.g. `mys1080-v4`) takes
an `[outprefix]` before those:

    ./dym-mod-metro-savecfg <group> D Nt Nx beta0 beta1 beta2 seed [outprefix] [K] [N] [Ntherm]

## Analysis

The Python scripts in `scripts/` parse the `GMES:` lines and produce plaquette
("freeze transition") curves, `(beta1, beta2)` contour/phase plots, and
Polyakov-loop diagnostics.

> Note: several analysis scripts were written with hard-coded absolute data paths
> (e.g. `/home/guest/dym_par_adj/...`). Update those to point at this repo's
> `data/` folder before rerunning.
