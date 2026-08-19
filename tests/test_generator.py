import numpy as np
import pytest

from ras_eif.generator import PopulationTransfer, build_generator, validate_generator


def test_generator_conserves_population_by_column() -> None:
    generator = build_generator(
        3,
        [
            PopulationTransfer(0, 1, 2.0),
            PopulationTransfer(1, 0, 1.0),
            PopulationTransfer(1, 2, 3.0),
        ],
    )
    assert np.allclose(generator.sum(axis=0), 0.0)
    assert generator[1, 0] == 2.0
    assert generator[0, 0] == -2.0
    validate_generator(generator)


def test_generator_rejects_negative_off_diagonal_rate() -> None:
    generator = np.array([[-1.0, -0.1], [1.0, 0.1]])
    with pytest.raises(ValueError, match="negative off-diagonal"):
        validate_generator(generator)


def test_generator_conservation_check_scales_with_rate_magnitude() -> None:
    generator = build_generator(
        2,
        [PopulationTransfer(0, 1, 1.0e6), PopulationTransfer(1, 0, 2.0e6)],
    )
    validate_generator(generator)
