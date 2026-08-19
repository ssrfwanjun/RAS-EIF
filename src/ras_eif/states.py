"""Atomic-state identifiers independent of any array ordering."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

_L_LABELS = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}


def parse_l(value: int | str) -> int:
    """Return an orbital quantum number from an integer or spectroscopic label."""
    if isinstance(value, int):
        return value
    try:
        return _L_LABELS[value.upper()]
    except KeyError as error:
        raise ValueError(f"Unsupported orbital label: {value!r}") from error


def parse_j(value: float | str | Fraction) -> Fraction:
    """Return an exact total-angular-momentum quantum number."""
    return Fraction(str(value))


@dataclass(frozen=True, order=True)
class AtomicState:
    """Fine-structure atomic state identified by physical quantum numbers."""

    n: int
    l: int
    j: Fraction

    def __post_init__(self) -> None:
        if self.n < 1:
            raise ValueError("n must be positive")
        if not 0 <= self.l < self.n:
            raise ValueError(f"l={self.l} is invalid for n={self.n}")
        valid_j = {Fraction(1, 2)} if self.l == 0 else {Fraction(2 * self.l - 1, 2), Fraction(2 * self.l + 1, 2)}
        if self.j not in valid_j:
            raise ValueError(f"j={self.j} is invalid for n={self.n}, l={self.l}")

    @classmethod
    def from_config(cls, value: list[object] | tuple[object, object, object]) -> AtomicState:
        if len(value) != 3:
            raise ValueError("A state must be [n, l, j]")
        return cls(n=int(value[0]), l=parse_l(value[1]), j=parse_j(value[2]))

    @property
    def label(self) -> str:
        letters = "SPDFGHI"
        return f"{self.n}{letters[self.l]}{self.j.numerator}/{self.j.denominator}"
