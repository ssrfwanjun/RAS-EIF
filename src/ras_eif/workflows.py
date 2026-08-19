"""High-level, explicitly scoped workflows for EIF visualizations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .arc_backend import ArcCesiumBackend
from .config import TransitionConfig, simulation_for_plot, states_for_simulation
from .generator import PopulationTransfer, build_generator
from .model import (
    atomic_frequencies_hz,
    driven_pair_transfers,
    ground_connected_indices,
    radiative_transfers,
    resolve_state_indices,
    spontaneous_rate_matrix,
)
from .spectrum import EmissionLine, bin_lines, emission_lines
from .steady_state import solve_steady_state


@dataclass(frozen=True)
class EifSimulation:
    """Results from a declared qualitative EIF rate-model calculation."""

    populations: np.ndarray
    frequencies_hz: np.ndarray
    spontaneous_rates_s: np.ndarray
    lines: tuple[EmissionLine, ...]


def simulate_qualitative_eif(
    config: TransitionConfig, backend: ArcCesiumBackend, thz_on: bool
) -> EifSimulation:
    """Calculate a qualitative EIF spectrum for one configured transition.

    Any spontaneous decay rate not represented by a state in `state_space` is
    returned to the ground-state reservoir. This
    prevents hidden population loss but also means the result is a visualization,
    not a converged all-channel calculation.
    """
    simulation = simulation_for_plot(config)
    candidate_states = list(states_for_simulation(config))
    candidate_frequencies = atomic_frequencies_hz(candidate_states, backend)
    # ARC exposes some low-n numerical states below the Cs 6S ground state.
    # They cannot be population sinks in this valence-electron model.
    states = [
        state for state, frequency in zip(candidate_states, candidate_frequencies) if frequency >= -1.0
    ]
    frequencies = candidate_frequencies[candidate_frequencies >= -1.0]
    if not set(config.ladder).issubset(states):
        raise ValueError("The configured excitation ladder contains a state below the Cs 6S ground state")
    indices = resolve_state_indices(states)
    spontaneous = spontaneous_rate_matrix(states, backend, frequencies)
    retained = ground_connected_indices(spontaneous, indices[config.ladder[0]])
    states = [states[index] for index in retained]
    frequencies = frequencies[retained]
    spontaneous = spontaneous[np.ix_(retained, retained)]
    if not set(config.ladder).issubset(states):
        raise ValueError("The configured ladder is not radiatively connected to the Cs 6S ground state")
    indices = resolve_state_indices(states)
    transfers = radiative_transfers(spontaneous)
    ground = config.ladder[0]
    ground_index = indices[ground]

    for upper, state in enumerate(states):
        if upper == ground_index:
            continue
        total_rate = 1.0 / backend.spontaneous_lifetime_s(state)
        represented_rate = spontaneous[upper, :].sum()
        unresolved_rate = max(total_rate - represented_rate, 0.0)
        if unresolved_rate:
            transfers.append(
                PopulationTransfer(upper, ground_index, unresolved_rate, "unresolved radiative reservoir")
            )

    for drive in simulation.drives:
        if drive.thz_only and not thz_on:
            continue
        if drive.effective_rate_s is not None:
            lower = indices[drive.lower]
            upper = indices[drive.upper]
            transfers.extend(
                (
                    PopulationTransfer(lower, upper, drive.effective_rate_s, "effective driven absorption"),
                    PopulationTransfer(upper, lower, drive.effective_rate_s, "effective driven emission"),
                )
            )
        else:
            transfers.extend(
                driven_pair_transfers(
                    states,
                    drive.lower,
                    drive.upper,
                    backend,
                    drive.intensity_w_m2 or 0.0,
                    drive.dephasing_rad_s or 0.0,
                    drive.detuning_rad_s,
                    drive.angular_factor,
                )
            )

    generator = build_generator(len(states), transfers)
    populations = solve_steady_state(generator)
    return EifSimulation(
        populations=populations,
        frequencies_hz=frequencies,
        spontaneous_rates_s=spontaneous,
        lines=tuple(emission_lines(populations, frequencies, spontaneous)),
    )


def plot_qualitative_eif(
    config: TransitionConfig,
    thz_off: EifSimulation,
    thz_on: EifSimulation,
    output: str | Path,
) -> Path:
    """Save a source-code-style THz-off/on/EIF stick-spectrum figure."""
    simulation = simulation_for_plot(config)
    wavelength = np.linspace(*simulation.wavelength_window_nm, 20_001)
    off = bin_lines(thz_off.lines, wavelength)
    on = bin_lines(thz_on.lines, wavelength)
    induced = on - off
    scale = max(np.max(np.abs(off)), np.max(np.abs(on)), np.max(np.abs(induced)))
    if scale == 0:
        raise RuntimeError("No emission lines fall within the configured wavelength window")
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8), constrained_layout=True)
    panels = ((off, "THz off", "#c33"), (on, "THz on", "#222"), (induced, "THz-induced EIF (on - off)", "#147d79"))
    for axis, (spectrum, title, color) in zip(axes, panels):
        axis.vlines(wavelength, 0.0, spectrum / scale, colors=color, linewidth=0.8)
        axis.axhline(0.0, color="#555", linewidth=0.5)
        axis.set_ylabel("Relative power")
        axis.set_title(title, loc="left", fontsize=11)
    axes[-1].set_xlabel("Wavelength (nm)")
    figure.suptitle(f"{config.identifier}: qualitative rate-model visualization", fontsize=13)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)
    return output_path
