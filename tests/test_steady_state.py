import numpy as np

from ras_eif.generator import PopulationTransfer, build_generator
from ras_eif.steady_state import solve_steady_state


def test_steady_state_replaces_redundant_equation_with_normalization() -> None:
    generator = build_generator(
        2,
        [
            PopulationTransfer(0, 1, 2.0),
            PopulationTransfer(1, 0, 1.0),
        ],
    )
    population = solve_steady_state(generator, total_population=12.0)
    assert np.allclose(generator @ population, 0.0)
    assert np.isclose(population.sum(), 12.0)
    assert np.allclose(population, [4.0, 8.0])


def test_raw_generator_is_singular_but_augmented_solution_is_not() -> None:
    generator = build_generator(2, [PopulationTransfer(0, 1, 1.0)])
    assert np.linalg.matrix_rank(generator) == 1
    population = solve_steady_state(generator)
    assert np.allclose(population, [0.0, 1.0])


def test_steady_state_accepts_scaled_rates_with_a_relative_residual_check() -> None:
    generator = build_generator(
        2,
        [
            PopulationTransfer(0, 1, 1.0e6),
            PopulationTransfer(1, 0, 2.0e6),
        ],
    )
    population = solve_steady_state(generator)
    assert np.allclose(population, [2.0 / 3.0, 1.0 / 3.0])
