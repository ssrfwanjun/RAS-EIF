from fractions import Fraction

import numpy as np
import pytest

from ras_eif.model import ground_connected_indices, resolve_state_indices
from ras_eif.rates import effective_incoherent_rate
from ras_eif.states import AtomicState


def test_state_index_map_rejects_duplicate_physical_states() -> None:
    ground = AtomicState(6, 0, Fraction(1, 2))
    with pytest.raises(ValueError, match="duplicate"):
        resolve_state_indices([ground, ground])


def test_effective_rate_decreases_away_from_resonance() -> None:
    on_resonance = effective_incoherent_rate(10.0, 5.0, 0.0)
    detuned = effective_incoherent_rate(10.0, 5.0, 20.0)
    assert on_resonance > detuned > 0.0


def test_ground_connected_indices_excludes_radiative_dead_ends() -> None:
    # 2 -> 1 -> 0 reaches ground; state 3 is an artificial terminal state.
    rates = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )
    assert np.array_equal(ground_connected_indices(rates, 0), [0, 1, 2])
