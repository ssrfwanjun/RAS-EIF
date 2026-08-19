import numpy as np

from ras_eif.spectrum import (
    PLANCK_CONSTANT,
    SPEED_OF_LIGHT,
    EmissionLine,
    bin_lines,
    convolve_gaussian,
    emission_lines,
)


def test_emission_line_power_is_population_times_a_h_nu() -> None:
    populations = np.array([0.7, 0.3])
    frequencies = np.array([0.0, 6.0e14])
    rates = np.array([[0.0, 0.0], [2.0e6, 0.0]])
    lines = emission_lines(populations, frequencies, rates)
    assert len(lines) == 1
    assert np.isclose(lines[0].power_w, 0.3 * 2.0e6 * PLANCK_CONSTANT * 6.0e14)
    assert np.isclose(lines[0].wavelength_nm, SPEED_OF_LIGHT / 6.0e14 * 1e9)


def test_gaussian_convolution_preserves_integrated_line_power() -> None:
    populations = np.array([0.0, 1.0])
    frequencies = np.array([0.0, 6.0e14])
    rates = np.array([[0.0, 0.0], [1.0, 0.0]])
    lines = emission_lines(populations, frequencies, rates)
    grid = np.linspace(lines[0].wavelength_nm - 2.0, lines[0].wavelength_nm + 2.0, 10001)
    spectrum = convolve_gaussian(lines, grid, fwhm_nm=0.04)
    assert np.isclose(np.trapezoid(spectrum, grid), lines[0].power_w, rtol=1e-6)


def test_binning_discards_lines_outside_the_display_window() -> None:
    wavelength = np.linspace(450.0, 650.0, 11)
    lines = (
        EmissionLine(1, 0, 500.0, 2.0, 1.0),
        EmissionLine(1, 0, 700.0, 9.0, 1.0),
    )
    spectrum = bin_lines(lines, wavelength)
    assert np.isclose(spectrum.sum(), 2.0)

