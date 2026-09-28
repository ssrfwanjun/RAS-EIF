# Legacy-vs-effective rate-model example

This example compares the legacy relative-rate heuristic with the revised
intensity-based effective-rate calculation for `24D32_26P12`.

## Scope

Both models use the revised named-state cascade manifold, ARC 3.10.2
zero-temperature spontaneous rates, the same unresolved-radiative reservoir,
and the same steady-state solver. Only the driven rates are changed. The
legacy formula is reconstructed using named states and ARC reduced dipoles; the
archived positional-index script is not executed.

The parameters are illustrative and are not experimental calibration values.
The effective-rate baseline is calibrated only so that its first drive equals
the legacy reference rate of `1.0e6 s^-1`; this makes the comparison isolate
the effects of the angular factor, dephasing, and detuning choices.

## Parameters

```yaml
transition: 24D32_26P12
state_space: n=4..26 cascade
legacy_reference_rate_s: 1000000.0
legacy_powers_mw: [10.0, 10.0, 1.0]
dephasing_rad_s: 100000000.0
detuning_over_dephasing: [0.0, 0.05, 0.05]
angular_factor: [1.0, 0.98, 0.98]
```

## Drive-rate comparison

| Drive | Transition | Legacy off (s⁻¹) | Effective off (s⁻¹) | Legacy on (s⁻¹) | Effective on (s⁻¹) | On-rate change |
|---|---|---:|---:|---:|---:|---:|
| W1 | 6S1/2 → 6P3/2 | 1e+06 | 1e+06 | 1e+06 | 1e+06 | +0.000% |
| W2 | 6P3/2 → 24D3/2 | 88.5016 | 84.1554 | 88.5016 | 84.1554 | -4.911% |
| W3 | 24D3/2 → 26P1/2 | 0 | 0 | 1.60723e+08 | 1.5283e+08 | -4.911% |

## Spectrum comparison

| Spectrum | Max difference / common max | L1 difference / legacy L1 |
|---|---:|---:|
| off | 4.9081% | 4.9081% |
| on | 4.9083% | 4.9083% |
| induced | 4.9067% | 4.9076% |

The spectra are convolved with the configured 0.04 nm Gaussian instrument
response. The corresponding spectral arrays are in
`spectrum_comparison.csv`; the plot uses a common normalization across both
models and is in `rate_model_comparison.png`.

## Population comparison

| Condition | L1 population difference | Maximum absolute difference |
|---|---:|---:|
| off | 1.85713e-06 | 8.99206e-07 |
| on | 3.07297e-06 | 1.48967e-06 |

## Interpretation

The numerical differences reported above quantify this matched illustrative
example only. They support a statement about the stability of the qualitative
rate-model visualization under this particular refinement; they do not replace
transition-specific calibration, polarization-resolved modeling, or state-space
convergence for a quantitative publication calculation.

## Reproduction

```bash
uv run python scripts/compare_rate_models.py --output-dir result/example
```
