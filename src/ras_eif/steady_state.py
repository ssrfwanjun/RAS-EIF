"""Steady-state solvers for conservative population-rate equations."""

from __future__ import annotations

import numpy as np

from .generator import validate_generator


def solve_steady_state(
    generator: np.ndarray,
    total_population: float = 1.0,
    normalization_row: int = -1,
    atol: float = 1e-10,
) -> np.ndarray:
    """Solve `R p = 0` subject to `sum(p) = total_population`.

    The raw generator is deliberately singular. A copied matrix has one
    redundant row replaced by the normalization constraint, making the system
    solvable without inverting the raw generator.
    """
    if total_population <= 0:
        raise ValueError("total_population must be positive")
    validate_generator(generator, atol=atol)
    state_count = generator.shape[0]
    row = normalization_row % state_count
    augmented = generator.copy()
    rhs = np.zeros(state_count, dtype=float)
    augmented[row, :] = 1.0
    rhs[row] = total_population
    population = np.linalg.solve(augmented, rhs)
    residual = np.linalg.norm(generator @ population, ord=np.inf)
    scale = max(np.linalg.norm(generator, ord=np.inf) * np.linalg.norm(population, ord=np.inf), 1.0)
    if residual > atol * scale:
        raise RuntimeError("Steady-state residual exceeds tolerance")
    if not np.isclose(population.sum(), total_population, atol=atol):
        raise RuntimeError("Steady-state normalization failed")
    if np.any(population < -atol):
        raise RuntimeError("Steady-state solution has negative populations")
    return np.maximum(population, 0.0)
