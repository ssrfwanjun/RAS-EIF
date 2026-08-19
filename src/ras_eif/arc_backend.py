"""ARC adapter that resolves states by quantum numbers, never list position."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .states import AtomicState


@dataclass
class ArcCesiumBackend:
    """Cesium data access with ARC imported only when a calculation is requested."""

    temperature_k: float = 0.0

    def __post_init__(self) -> None:
        try:
            from arc import Cesium
        except ImportError as error:  # pragma: no cover - depends on local ARC installation
            raise RuntimeError(
                "ARC is required for atomic-data calculations. Install "
                "ARC-Alkali-Rydberg-Calculator."
            ) from error
        self.atom = Cesium()

    @staticmethod
    def _args(state: AtomicState) -> tuple[int, int, float]:
        return state.n, state.l, float(state.j)

    def transition_frequency_hz(self, lower: AtomicState, upper: AtomicState) -> float:
        """Return signed transition frequency calculated by ARC."""
        lower_args = self._args(lower)
        upper_args = self._args(upper)
        return self.atom.getTransitionFrequency(*lower_args, *upper_args)

    def frequency_from_ground_hz(self, state: AtomicState) -> float:
        ground = AtomicState(6, 0, Fraction(1, 2))
        return self.transition_frequency_hz(ground, state)

    def reduced_dipole_au(self, first: AtomicState, second: AtomicState) -> float:
        """Return the magnitude of ARC's reduced J-basis dipole matrix element."""
        return abs(self.atom.getReducedMatrixElementJ(*self._args(first), *self._args(second)))

    def spontaneous_rate_s(self, upper: AtomicState, lower: AtomicState) -> float:
        """Return the zero-temperature downward rate, or zero for non-E1 channels."""
        if upper.l - lower.l not in (-1, 1):
            return 0.0
        frequency = self.transition_frequency_hz(lower, upper)
        if frequency <= 0:
            return 0.0
        return self.atom.getTransitionRate(*self._args(upper), *self._args(lower), temperature=0.0)

    def spontaneous_lifetime_s(self, state: AtomicState) -> float:
        """Return ARC's zero-temperature state lifetime in seconds."""
        lifetime = self.atom.getStateLifetime(*self._args(state), temperature=0.0)
        if lifetime <= 0:
            raise ValueError(f"ARC returned a non-positive lifetime for {state.label}")
        return lifetime


def generate_fine_structure_states(n_min: int, n_max: int) -> list[AtomicState]:
    """Generate all fine-structure states in the declared principal-number range."""
    if n_min < 1 or n_max <= n_min:
        raise ValueError("Require 1 <= n_min < n_max")
    states: list[AtomicState] = []
    for n in range(n_min, n_max):
        for l in range(n):
            if l == 0:
                states.append(AtomicState(n, l, Fraction(1, 2)))
            else:
                states.extend((AtomicState(n, l, Fraction(2 * l - 1, 2)), AtomicState(n, l, Fraction(2 * l + 1, 2))))
    return states
