"""Assembly helpers that connect atomic data to the validated rate solver."""

from __future__ import annotations

import numpy as np

from .arc_backend import ArcCesiumBackend
from .generator import PopulationTransfer
from .rates import effective_incoherent_rate, rabi_frequency
from .states import AtomicState


def resolve_state_indices(states: list[AtomicState]) -> dict[AtomicState, int]:
    """Return a unique state-to-index map and reject duplicate physical states."""
    indices = {state: index for index, state in enumerate(states)}
    if len(indices) != len(states):
        raise ValueError("State space contains duplicate physical states")
    return indices


def atomic_frequencies_hz(states: list[AtomicState], backend: ArcCesiumBackend) -> np.ndarray:
    """Calculate frequencies relative to 6S1/2 for every declared state."""
    return np.array([backend.frequency_from_ground_hz(state) for state in states], dtype=float)


def spontaneous_rate_matrix(
    states: list[AtomicState], backend: ArcCesiumBackend, frequencies_hz: np.ndarray | None = None
) -> np.ndarray:
    """Return A[upper, lower] for all allowed downward channels in the state space.

    Channels outside the declared state space are deliberately not hidden. A
    production calculation must either expand the state space until convergence
    or add an explicit, documented loss sink.
    """
    frequencies = frequencies_hz if frequencies_hz is not None else atomic_frequencies_hz(states, backend)
    if frequencies.shape != (len(states),):
        raise ValueError("frequencies_hz has the wrong shape")
    rates = np.zeros((len(states), len(states)), dtype=float)
    for upper, upper_state in enumerate(states):
        for lower, lower_state in enumerate(states):
            if frequencies[upper] > frequencies[lower]:
                rates[upper, lower] = backend.spontaneous_rate_s(upper_state, lower_state)
    return rates


def radiative_transfers(spontaneous_rates_s: np.ndarray) -> list[PopulationTransfer]:
    """Convert A[upper, lower] into population transfers for a rate generator."""
    count = spontaneous_rates_s.shape[0]
    if spontaneous_rates_s.shape != (count, count):
        raise ValueError("spontaneous_rates_s must be square")
    return [
        PopulationTransfer(source=upper, destination=lower, rate_s=float(rate), label="spontaneous")
        for upper in range(count)
        for lower in range(count)
        if (rate := spontaneous_rates_s[upper, lower]) > 0.0
    ]


def ground_connected_indices(spontaneous_rates_s: np.ndarray, ground_index: int) -> np.ndarray:
    """Return states that can reach the ground through downward radiative channels."""
    count = spontaneous_rates_s.shape[0]
    if spontaneous_rates_s.shape != (count, count):
        raise ValueError("spontaneous_rates_s must be square")
    if not 0 <= ground_index < count:
        raise IndexError("ground_index is outside the state space")
    reachable = np.zeros(count, dtype=bool)
    reachable[ground_index] = True
    while True:
        # A state can reach the ground if it decays to a state already known to
        # reach the ground. Rows are upper states and columns are lower states.
        newly_reachable = np.any(spontaneous_rates_s[:, reachable] > 0.0, axis=1)
        updated = reachable | newly_reachable
        if np.array_equal(updated, reachable):
            return np.flatnonzero(reachable)
        reachable = updated


def driven_pair_transfers(
    states: list[AtomicState],
    lower: AtomicState,
    upper: AtomicState,
    backend: ArcCesiumBackend,
    intensity_w_m2: float,
    dephasing_rad_s: float,
    detuning_rad_s: float = 0.0,
    angular_factor: float = 1.0,
) -> list[PopulationTransfer]:
    """Build bidirectional effective transfers for one declared driven transition."""
    indices = resolve_state_indices(states)
    if lower not in indices or upper not in indices:
        raise KeyError("Driven states must both be included in the state space")
    dipole = backend.reduced_dipole_au(lower, upper)
    omega = rabi_frequency(intensity_w_m2, dipole, angular_factor)
    rate = effective_incoherent_rate(omega, dephasing_rad_s, detuning_rad_s)
    return [
        PopulationTransfer(indices[lower], indices[upper], rate, "driven absorption"),
        PopulationTransfer(indices[upper], indices[lower], rate, "driven stimulated emission"),
    ]
