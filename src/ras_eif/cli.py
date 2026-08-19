"""Command-line entry points for configuration validation."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from .arc_backend import ArcCesiumBackend
from .config import load_transition_config, simulation_for_plot
from .workflows import plot_qualitative_eif, simulate_qualitative_eif


def _validate_transition(path: str, tolerance_mhz: float) -> int:
    config = load_transition_config(path)
    backend = ArcCesiumBackend()
    print(f"Transition: {config.identifier}")
    for lower, upper in zip(config.ladder, config.ladder[1:]):
        frequency_ghz = backend.transition_frequency_hz(lower, upper) / 1e9
        print(f"  {lower.label} -> {upper.label}: {frequency_ghz:.9f} GHz")
    calculated_ghz = backend.transition_frequency_hz(config.ladder[-2], config.ladder[-1]) / 1e9
    difference_mhz = abs(calculated_ghz - config.expected_thz_ghz) * 1e3
    print(f"Expected THz: {config.expected_thz_ghz:.9f} GHz")
    print(f"Difference: {difference_mhz:.3f} MHz")
    if difference_mhz > tolerance_mhz:
        print(f"FAILED: difference exceeds {tolerance_mhz:.3f} MHz")
        return 1
    print("PASSED: state identities and final transition frequency are consistent.")
    return 0


def _plot_spectrum(path: str, output: str | None) -> int:
    config = load_transition_config(path)
    if config.simulation is None:
        print("Using default qualitative simulation inputs; no calibrated simulation section was supplied.")
    backend = ArcCesiumBackend()
    thz_off = simulate_qualitative_eif(config, backend, thz_on=False)
    thz_on = simulate_qualitative_eif(config, backend, thz_on=True)
    lower_wavelength, upper_wavelength = simulation_for_plot(config).wavelength_window_nm
    off_count = sum(lower_wavelength <= line.wavelength_nm <= upper_wavelength for line in thz_off.lines)
    on_count = sum(lower_wavelength <= line.wavelength_nm <= upper_wavelength for line in thz_on.lines)
    destination = Path(output) if output else Path("result/generated") / f"{config.identifier}_eif.png"
    saved = plot_qualitative_eif(config, thz_off, thz_on, destination)
    print(f"Emission channels in {lower_wavelength:g}-{upper_wavelength:g} nm: THz off={off_count}, THz on={on_count}")
    print(f"Saved qualitative EIF spectrum: {saved}")
    return 0


def _validate_directory(directory: str, tolerance_mhz: float) -> int:
    paths = sorted(Path(directory).rglob("*.yaml"))
    if not paths:
        raise ValueError(f"No YAML files found in {directory}")
    failures = 0
    for path in paths:
        print(f"\n== {path} ==")
        failures += _validate_transition(str(path), tolerance_mhz)
    print(f"\nValidated {len(paths)} transition configurations; failures: {failures}")
    return int(failures > 0)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ras-eif")
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate", help="validate an explicitly named transition")
    validate.add_argument("config", help="path to a transition YAML file")
    validate.add_argument("--tolerance-mhz", type=float, default=5.0)
    validate_all = subparsers.add_parser("validate-all", help="validate every transition YAML in a directory")
    validate_all.add_argument("directory", help="directory containing transition YAML files")
    validate_all.add_argument("--tolerance-mhz", type=float, default=5.0)
    spectrum = subparsers.add_parser("spectrum", help="plot configured qualitative EIF stick spectra")
    spectrum.add_argument("config", help="path to a transition YAML file")
    spectrum.add_argument("--output", help="PNG output path")
    arguments = parser.parse_args(argv)
    if arguments.command == "validate":
        try:
            return _validate_transition(arguments.config, arguments.tolerance_mhz)
        except RuntimeError as error:
            parser.exit(2, f"ERROR: {error}\n")
    if arguments.command == "spectrum":
        try:
            return _plot_spectrum(arguments.config, arguments.output)
        except (RuntimeError, ValueError) as error:
            parser.exit(2, f"ERROR: {error}\n")
    if arguments.command == "validate-all":
        try:
            return _validate_directory(arguments.directory, arguments.tolerance_mhz)
        except (RuntimeError, ValueError) as error:
            parser.exit(2, f"ERROR: {error}\n")
    parser.error("Unknown command")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
