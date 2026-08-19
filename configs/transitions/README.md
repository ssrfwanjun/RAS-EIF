# Transition Configurations

Each YAML file replaces one legacy Python script and identifies the excitation
ladder by physical `(n, l, j)` quantum numbers. The file name and `id` match the
legacy script stem so that the migration is auditable.

The directory groups all 31 archived channels by excitation scheme:

- `two-photon/`: 19 paths, `6S1/2 -> 6P3/2 -> nS/nD -> nP`.
- `three-photon/`: 12 paths, `6S1/2 -> 6P3/2 -> 7S1/2 -> nP -> nD`.

Validate every path against the locked ARC release with:

```bash
uv run ras-eif validate-all configs/transitions
```

`two-photon/24D32_26P12.yaml` additionally contains an explicitly marked qualitative
visualization section used by the `spectrum` command. Do not copy its effective
rates to another transition or use them for quantitative theory without
transition-specific calibration and state-space convergence tests.

All three-photon channels define an explicit four-drive `simulation` section.
Two-photon channels without an explicit section use the same full cascade
manifold from `n=4` through the highest ladder state, with default effective
rates until transition-specific drive conditions are supplied.

Each `simulation.drives` entry can use one of two modes:

- `effective_rate_s`: a declared phenomenological population-transfer rate.
- `intensity_w_m2` and `dephasing_rad_s`, optionally `detuning_rad_s` and
  `angular_factor`: an effective rate calculated from the ARC dipole element.

Do not specify both modes for one drive.

For a complete model-space calculation, configure the state manifold as:

```yaml
state_space:
  mode: cascade
  n_min: 4
  n_max: 47
```

This includes every fine-structure state in the stated inclusive principal
quantum-number range, in addition to the explicit ladder. It is substantially
slower than the reduced state space but is required to produce cascade lines.
States whose ARC energy lies below the Cs `6S1/2` ground state are excluded
before the rate generator is assembled.
