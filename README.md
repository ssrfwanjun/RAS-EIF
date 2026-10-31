# RAS-EIF

Reproducible rate-equation tools for modelling terahertz-induced
electromagnetically induced fluorescence (EIF) in cesium Rydberg atoms.

## Status

The scripts in `code/` are retained as an **unvalidated historical archive**.
They use position-based atomic-state indices and must not be used to generate
figures or support quantitative claims. New calculations belong in `src/`, are
configured by quantum numbers in `configs/`, and must pass the checks in
`tests/`.

This repository currently reconstructs the atomic rate-equation model only. It
does not yet reproduce the experimental acquisition, deep-learning classifier,
COMSOL simulation, or metasurface Fano fitting.

## Install

```bash
uv sync --extra dev
```

ARC 3.10.2 is a locked runtime dependency and is installed by this command.

## Validate a configuration

```bash
uv run ras-eif validate configs/transitions/two-photon/24D32_26P12.yaml
```

The command resolves each state by its explicit `(n, l, j)` quantum numbers and
reports the ladder transition frequencies. It never depends on a sorted-array
position.

Validate all archived excitation channels with:

```bash
uv run ras-eif validate-all configs/transitions
```

## Plot a qualitative EIF spectrum

```bash
uv run ras-eif spectrum configs/transitions/two-photon/24D32_26P12.yaml
```

This writes `result/generated/24D32_26P12_eif.png` with THz-off, THz-on, and
THz-induced stick spectra, in the same form as the legacy scripts. Its rates
are illustrative only; do not use this reduced-state figure as a publication
calculation before completing state-space convergence and parameter calibration.

All archived transition YAML files can be plotted. Files without a `simulation`
section use documented default effective rates and a full cascade manifold from
`n=4` through their highest ladder state; the CLI prints a warning whenever
those qualitative defaults are used.

## Generate all refactored simulation results

To calculate every archived channel with the refactored model in both THz-off
and THz-on conditions:

```bash
uv run python scripts/generate_new_simulation_results.py
```

This writes two MATLAB v5 files per channel under
`result/simulation/new/two-photon/` and `result/simulation/new/three-photon/`,
using the calibrated wavelength grid from `result/Frequency_Calibration.mat`.
The output arrays contain wavelength in nm and unit-total-population emitted
power in watts. The batch command uses ARC rates and state lifetimes at 300 K by
default and scales THz-only drives by 10; pass `--temperature-k 0` for a
zero-temperature run or `--thz-intensity-scale 1` for the unscaled THz run. The
run manifest records the temperature, THz scale, input mode, and population
normalization.

Compare the overlapping archived and refactored channels with normalized
shape metrics:

```bash
uv run python scripts/compare_old_new_simulation.py
```

The report is written to `result/simulation/new/comparison_old_vs_new.md` and
the row-level metrics to `comparison_old_vs_new.csv`.

Generate visual old/new overlays for every new channel with:

```bash
uv run python scripts/plot_old_new_simulation.py
```

This writes per-channel three-panel figures and condition-specific overview
figures under `result/simulation/new/plots/`.

## Model conventions

For a population column vector `p`, the raw rate generator is defined by

`dp/dt = R p`.

For each transition `source -> destination` with rate `k`, the generator adds
`R[destination, source] += k` and `R[source, source] -= k`. Thus every column
of `R` sums to zero. The steady state is found by replacing one redundant row
with `sum(p) = total_population`; the raw generator remains available and is
validated before that replacement.

The model uses fixed atomic Einstein coefficients. Applied fields redistribute
populations through effective transition rates; they do not change branching
ratios or Einstein coefficients.

## Layout

- `src/ras_eif/`: model implementation.
- `configs/transitions/`: explicitly named excitation ladders.
- `scripts/`: non-interactive validation and reproduction entry points.
- `tests/`: mathematical and configuration regression tests.
- `docs/MIGRATION.md`: known defects in the legacy scripts and migration rules.
- `result/`: retained source artifacts. Files under `result/experiment/` need
  provenance review before being described as experimental data.
- `code/`: legacy scripts.

## Before publishing regenerated spectra

1. Pin the ARC version and record all laser/THz intensities, beam areas,
   polarizations, linewidths, detunings, and temperatures.
2. Include all material radiative channels or quantify the omitted branching
   fraction using an explicit loss sink.
3. Convolve calculated line spectra with the measured instrument response.
4. Record the configuration, package versions, and validation diagnostics with
   every generated result.
