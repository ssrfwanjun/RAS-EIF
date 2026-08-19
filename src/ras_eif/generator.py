"""Construction and validation of conservative population-rate generators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PopulationTransfer:
    """One population transfer from `source` to `destination` at `rate_s`."""

    source: int
    destination: int
    rate_s: float
    label: str = ""

    def __post_init__(self) -> None:
        if self.source < 0 or self.destination < 0:
            raise ValueError("State indices must be non-negative")
        if self.source == self.destination:
            raise ValueError("A transfer must connect distinct states")
        if not np.isfinite(self.rate_s) or self.rate_s < 0:
            raise ValueError("Transfer rate must be finite and non-negative")


def build_generator(state_count: int, transfers: list[PopulationTransfer]) -> np.ndarray:
    """Build R for `dp/dt = R p` using source-column convention."""
    if state_count < 2:
        raise ValueError("At least two states are required")
    generator = np.zeros((state_count, state_count), dtype=float)
    for transfer in transfers:
        if transfer.source >= state_count or transfer.destination >= state_count:
            raise IndexError(f"Transfer {transfer.label!r} refers to a missing state")
        generator[transfer.destination, transfer.source] += transfer.rate_s
        generator[transfer.source, transfer.source] -= transfer.rate_s
    return generator


def validate_generator(generator: np.ndarray, atol: float = 1e-10) -> None:
    """Raise when a generator violates conservation or positivity conditions."""
    if generator.ndim != 2 or generator.shape[0] != generator.shape[1]:
        raise ValueError("Generator must be a square matrix")
    if not np.isfinite(generator).all():
        raise ValueError("Generator contains non-finite values")
    off_diagonal = generator.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if np.any(off_diagonal < -atol):
        raise ValueError("Generator has a negative off-diagonal rate")
    conservation_error = np.max(np.abs(generator.sum(axis=0)))
    scale = max(np.max(np.abs(generator)), 1.0)
    if conservation_error > atol * scale:
        raise ValueError("Generator does not conserve total population")
