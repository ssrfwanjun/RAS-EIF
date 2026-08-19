"""Rate and field-coupling utilities with explicit SI units."""

from __future__ import annotations

from math import sqrt

EPSILON_0 = 8.8541878128e-12
SPEED_OF_LIGHT = 299792458.0
HBAR = 1.054571817e-34
ELEMENTARY_CHARGE = 1.602176634e-19
BOHR_RADIUS = 5.29177210903e-11


def electric_field_from_intensity(intensity_w_m2: float) -> float:
    """Return the peak electric-field amplitude of a plane wave in V/m."""
    if intensity_w_m2 < 0:
        raise ValueError("Intensity must be non-negative")
    return sqrt(2.0 * intensity_w_m2 / (SPEED_OF_LIGHT * EPSILON_0))


def rabi_frequency(
    intensity_w_m2: float,
    dipole_au: float,
    angular_factor: float = 1.0,
) -> float:
    """Return |Omega| in rad/s from a reduced dipole and declared angular factor."""
    if not 0.0 <= angular_factor <= 1.0:
        raise ValueError("angular_factor must be between zero and one")
    dipole_si = abs(dipole_au) * ELEMENTARY_CHARGE * BOHR_RADIUS * angular_factor
    return dipole_si * electric_field_from_intensity(intensity_w_m2) / HBAR


def effective_incoherent_rate(
    omega_rad_s: float,
    dephasing_rad_s: float,
    detuning_rad_s: float = 0.0,
) -> float:
    """Return a documented population-transfer approximation in 1/s.

    This is an adiabatic rate-model approximation, not a replacement for an
    optical-Bloch calculation. The dephasing rate and detuning must therefore
    be experimentally justified before using this function quantitatively.
    """
    if dephasing_rad_s <= 0:
        raise ValueError("dephasing_rad_s must be positive")
    return omega_rad_s**2 * dephasing_rad_s / (dephasing_rad_s**2 + 4.0 * detuning_rad_s**2)


def thermal_photon_number(frequency_hz: float, temperature_k: float) -> float:
    """Return the blackbody photon occupation number."""
    if frequency_hz <= 0 or temperature_k <= 0:
        raise ValueError("frequency and temperature must be positive")
    exponent = 6.62607015e-34 * frequency_hz / (1.380649e-23 * temperature_k)
    if exponent > 700:
        return 0.0
    return 1.0 / (pow(2.718281828459045, exponent) - 1.0)
