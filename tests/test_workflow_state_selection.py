import numpy as np

from ras_eif.states import AtomicState


def test_ground_referenced_filter_excludes_states_below_6s() -> None:
    states = [
        AtomicState.from_config([6, "S", "1/2"]),
        AtomicState.from_config([5, "P", "3/2"]),
        AtomicState.from_config([7, "S", "1/2"]),
    ]
    frequencies = np.array([0.0, -1.0e14, 2.0e14])
    retained = [state for state, frequency in zip(states, frequencies) if frequency >= -1.0]
    assert [state.label for state in retained] == ["6S1/2", "7S1/2"]
