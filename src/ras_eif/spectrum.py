"""Emission-line generation and instrument-response convolution."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np

PLANCK_CONSTANT = 6.62607015e-34
SPEED_OF_LIGHT = 299792458.0


@dataclass(frozen=True)
class EmissionLine:
    """A spontaneous-emission line with power before collection losses."""

    upper: int
    lower: int
    wavelength_nm: float
    power_w: float
    rate_s: float


def emission_lines(
    populations: np.ndarray,
    frequencies_hz: np.ndarray,
    spontaneous_rates_s: np.ndarray,
) -> list[EmissionLine]:
    """Calculate `N_upper A_ul h nu` for supplied spontaneous channels only."""
    count = len(populations)
    if frequencies_hz.shape != (count,) or spontaneous_rates_s.shape != (count, count):
        raise ValueError("Population, frequency, and rate dimensions are inconsistent")
    lines: list[EmissionLine] = []
    for upper in range(count):
        for lower in range(count):
            rate = spontaneous_rates_s[upper, lower]
            frequency = frequencies_hz[upper] - frequencies_hz[lower]
            if rate <= 0.0 or frequency <= 0.0:
                continue
            wavelength_nm = SPEED_OF_LIGHT / frequency * 1e9
            power_w = populations[upper] * rate * PLANCK_CONSTANT * frequency
            lines.append(EmissionLine(upper, lower, wavelength_nm, power_w, rate))
    return lines


def convolve_gaussian(
    lines: list[EmissionLine], wavelength_nm: np.ndarray, fwhm_nm: float
) -> np.ndarray:
    """Return a spectrum convolved with a normalized Gaussian instrument response."""
    if fwhm_nm <= 0:
        raise ValueError("fwhm_nm must be positive")
    sigma = fwhm_nm / (2.0 * sqrt(2.0 * np.log(2.0)))
    spectrum = np.zeros_like(wavelength_nm, dtype=float)
    for line in lines:
        profile = np.exp(-0.5 * ((wavelength_nm - line.wavelength_nm) / sigma) ** 2)
        integral = np.trapezoid(profile, wavelength_nm)
        if integral > 0:
            spectrum += line.power_w * profile / integral
    return spectrum


def bin_lines(lines: list[EmissionLine] | tuple[EmissionLine, ...], wavelength_nm: np.ndarray) -> np.ndarray:
    """Bin line powers on a wavelength grid for source-code-style stick plots."""
    if wavelength_nm.ndim != 1 or len(wavelength_nm) < 2:
        raise ValueError("wavelength_nm must be a one-dimensional grid")
    spectrum = np.zeros_like(wavelength_nm, dtype=float)
    for line in lines:
        if not wavelength_nm[0] <= line.wavelength_nm <= wavelength_nm[-1]:
            continue
        index = np.searchsorted(wavelength_nm, line.wavelength_nm, side="right") - 1
        if 0 <= index < len(spectrum):
            spectrum[index] += line.power_w
    return spectrum

