# Legacy Code Migration

## Do not regenerate figures with `code/`

The files under `code/two-photon/` and `code/three-photon/` are retained only
to preserve the original deposited material. They contain duplicated scripts
and known correctness defects. In particular:

- They sort a list of states, delete rows, and then use stale integer indices
  for the excitation ladder.
- They explicitly invert an augmented matrix and subsequently normalize by the
  ground-state component instead of the total population.
- Their relative drive-rate scaling treats power as intensity and omits the
  beam area, polarization, linewidth, detuning, and angular factors.
- They use a mixture of 300 K transition rates and spontaneous-emission power
  without separating BBR-induced population transfer from fluorescence.
- They aggregate on/off spectral lines by floating-point wavelength equality.

## Replacement rules

1. Define every state as `(n, l, j)` in a YAML configuration.
2. Resolve states by quantum numbers and validate every ladder transition using
   ARC before calculation.
3. Build the raw conservative generator first. Test its off-diagonal positivity
   and zero column sum.
4. Solve the steady state only after replacing one redundant row with total
   population normalization.
5. Use zero-temperature spontaneous rates for emitted fluorescence. Add BBR,
   collisions, transit loss, and omitted channels as separate, documented
   transfers.
6. Store raw line lists and convolved spectra separately, with configuration,
   software version, units, and a declared `THz-on minus THz-off` convention.

## Scope of the new model

`src/ras_eif/` supplies a validated numerical foundation; it does not claim
that an incomplete state space is a complete physical model. Convergence with
respect to the included state space and treatment of collision-induced features
must be demonstrated for each publication figure.
