from fractions import Fraction

import pytest

from ras_eif.states import AtomicState, parse_l


def test_state_is_constructed_from_explicit_quantum_numbers() -> None:
    state = AtomicState.from_config([24, "D", "3/2"])
    assert state == AtomicState(24, 2, Fraction(3, 2))
    assert state.label == "24D3/2"


def test_invalid_quantum_numbers_are_rejected() -> None:
    with pytest.raises(ValueError):
        AtomicState.from_config([24, "D", "1/2"])
    with pytest.raises(ValueError):
        parse_l("Q")
