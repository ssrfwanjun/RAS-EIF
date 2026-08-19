"""Typed loading of explicit transition configurations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from .states import AtomicState


@dataclass(frozen=True)
class TransitionConfig:
    identifier: str
    ladder: tuple[AtomicState, ...]
    expected_thz_ghz: float
    figures: tuple[str, ...]
    simulation: SimulationConfig | None


@dataclass(frozen=True)
class EffectiveDrive:
    """One drive specified by an effective rate or physical field inputs."""

    lower: AtomicState
    upper: AtomicState
    effective_rate_s: float | None
    intensity_w_m2: float | None
    dephasing_rad_s: float | None
    detuning_rad_s: float
    angular_factor: float
    thz_only: bool


@dataclass(frozen=True)
class SimulationConfig:
    """Inputs for an auditable, qualitative EIF visualization."""

    state_space: tuple[AtomicState, ...]
    cascade_n_range: tuple[int, int] | None
    drives: tuple[EffectiveDrive, ...]
    wavelength_window_nm: tuple[float, float]
    instrument_fwhm_nm: float


def _load_simulation(raw: dict[object, object], ladder: tuple[AtomicState, ...]) -> SimulationConfig | None:
    simulation = raw.get("simulation")
    if simulation is None:
        return None
    if not isinstance(simulation, dict):
        raise TypeError("simulation must be a mapping")
    state_space = simulation["state_space"]
    cascade_n_range: tuple[int, int] | None = None
    if isinstance(state_space, dict):
        if state_space.get("mode") != "cascade":
            raise ValueError("state_space mapping mode must be 'cascade'")
        n_min = int(state_space["n_min"])
        n_max = int(state_space["n_max"])
        if n_min < 1 or n_max < n_min:
            raise ValueError("cascade state-space range is invalid")
        states = ladder
        cascade_n_range = (n_min, n_max)
    else:
        states = tuple(AtomicState.from_config(state) for state in state_space)
    if len(set(states)) != len(states):
        raise ValueError("simulation state_space contains duplicate states")
    if not set(ladder).issubset(states):
        raise ValueError("simulation state_space must contain the full excitation ladder")
    drives = tuple(_load_drive(item) for item in simulation["drives"])
    for drive in drives:
        if drive.lower not in states or drive.upper not in states:
            raise ValueError("Every driven state must be in simulation state_space")
    window = tuple(float(value) for value in simulation.get("wavelength_window_nm", [450, 650]))
    if len(window) != 2 or window[0] >= window[1]:
        raise ValueError("wavelength_window_nm must be [minimum, maximum]")
    return SimulationConfig(
        state_space=states,
        cascade_n_range=cascade_n_range,
        drives=drives,
        wavelength_window_nm=(window[0], window[1]),
        instrument_fwhm_nm=float(simulation.get("instrument_fwhm_nm", 0.04)),
    )


def _load_drive(raw: dict[object, object]) -> EffectiveDrive:
    rate = raw.get("effective_rate_s")
    intensity = raw.get("intensity_w_m2")
    if rate is not None and intensity is not None:
        raise ValueError("A drive must specify either effective_rate_s or intensity_w_m2, not both")
    if rate is None and intensity is None:
        raise ValueError("A drive must specify effective_rate_s or intensity_w_m2")
    dephasing = raw.get("dephasing_rad_s")
    if intensity is not None and dephasing is None:
        raise ValueError("An intensity-based drive requires dephasing_rad_s")
    drive = EffectiveDrive(
        lower=AtomicState.from_config(raw["lower"]),
        upper=AtomicState.from_config(raw["upper"]),
        effective_rate_s=float(rate) if rate is not None else None,
        intensity_w_m2=float(intensity) if intensity is not None else None,
        dephasing_rad_s=float(dephasing) if dephasing is not None else None,
        detuning_rad_s=float(raw.get("detuning_rad_s", 0.0)),
        angular_factor=float(raw.get("angular_factor", 1.0)),
        thz_only=bool(raw.get("thz_only", False)),
    )
    if drive.effective_rate_s is not None and drive.effective_rate_s <= 0:
        raise ValueError("effective_rate_s must be positive")
    if drive.intensity_w_m2 is not None and drive.intensity_w_m2 < 0:
        raise ValueError("intensity_w_m2 must be non-negative")
    if drive.dephasing_rad_s is not None and drive.dephasing_rad_s <= 0:
        raise ValueError("dephasing_rad_s must be positive")
    if not 0.0 <= drive.angular_factor <= 1.0:
        raise ValueError("angular_factor must be between zero and one")
    return drive


def load_transition_config(path: str | Path) -> TransitionConfig:
    """Load one transition configuration and reject ambiguous or incomplete ladders."""
    with Path(path).open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    ladder = tuple(AtomicState.from_config(state) for state in raw["ladder"])
    if len(ladder) < 3:
        raise ValueError("A ladder must contain ground, optically excited, and Rydberg states")
    if len(set(ladder)) != len(ladder):
        raise ValueError("A ladder cannot repeat a state")
    return TransitionConfig(
        identifier=str(raw["id"]),
        ladder=ladder,
        expected_thz_ghz=float(raw["expected_thz_ghz"]),
        figures=tuple(raw.get("figures", [])),
        simulation=_load_simulation(raw, ladder),
    )


def simulation_for_plot(config: TransitionConfig) -> SimulationConfig:
    """Return explicit simulation inputs or conservative qualitative defaults.

    Defaults make every archived channel plottable while keeping the YAML files
    focused on physical state identities. They are deliberately effective-rate
    visualizations, not measured or fitted experimental parameters.
    """
    if config.simulation is not None:
        return config.simulation
    drives = tuple(
        EffectiveDrive(
            lower=lower,
            upper=upper,
            effective_rate_s=1.0e6,
            intensity_w_m2=None,
            dephasing_rad_s=None,
            detuning_rad_s=0.0,
            angular_factor=1.0,
            thz_only=index == len(config.ladder) - 2,
        )
        for index, (lower, upper) in enumerate(zip(config.ladder, config.ladder[1:]))
    )
    return SimulationConfig(
        state_space=config.ladder,
        cascade_n_range=(4, max(state.n for state in config.ladder)),
        drives=drives,
        wavelength_window_nm=(450.0, 650.0),
        instrument_fwhm_nm=0.04,
    )


def states_for_simulation(config: TransitionConfig) -> tuple[AtomicState, ...]:
    """Expand a configured cascade manifold while preserving explicit ladder states."""
    simulation = simulation_for_plot(config)
    if simulation.cascade_n_range is None:
        return simulation.state_space
    from .arc_backend import generate_fine_structure_states

    n_min, n_max = simulation.cascade_n_range
    manifold = generate_fine_structure_states(n_min, n_max + 1)
    return tuple(dict.fromkeys((*config.ladder, *manifold, *simulation.state_space)))
