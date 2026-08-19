from pathlib import Path

from ras_eif.config import load_transition_config, simulation_for_plot, states_for_simulation


def test_transition_config_uses_named_quantum_states() -> None:
    root = Path(__file__).resolve().parents[1]
    config = load_transition_config(root / "configs/transitions/two-photon/24D32_26P12.yaml")
    assert config.identifier == "24D32_26P12"
    assert [state.label for state in config.ladder] == ["6S1/2", "6P3/2", "24D3/2", "26P1/2"]
    assert config.expected_thz_ghz == 548.700


def test_config_without_simulation_receives_qualitative_plot_defaults() -> None:
    root = Path(__file__).resolve().parents[1]
    config = load_transition_config(root / "configs/transitions/two-photon/29S12_30P12.yaml")
    simulation = simulation_for_plot(config)
    states = states_for_simulation(config)
    assert len(simulation.drives) == len(config.ladder) - 1
    assert simulation.drives[-1].thz_only
    assert simulation.state_space == config.ladder
    assert len(states) > 100


def test_intensity_drive_configuration_requires_dephasing_and_is_loaded(tmp_path: Path) -> None:
    path = tmp_path / "intensity.yaml"
    path.write_text(
        """
id: intensity_test
ladder: [[6, S, 1/2], [6, P, 3/2], [24, D, 3/2], [26, P, 1/2]]
expected_thz_ghz: 548.7
simulation:
  state_space: [[6, S, 1/2], [6, P, 3/2], [24, D, 3/2], [26, P, 1/2]]
  drives:
    - lower: [6, S, 1/2]
      upper: [6, P, 3/2]
      intensity_w_m2: 100.0
      dephasing_rad_s: 1.0e7
    - lower: [6, P, 3/2]
      upper: [24, D, 3/2]
      effective_rate_s: 1.0e6
    - lower: [24, D, 3/2]
      upper: [26, P, 1/2]
      effective_rate_s: 1.0e6
      thz_only: true
""",
        encoding="utf-8",
    )
    config = load_transition_config(path)
    drive = config.simulation.drives[0] if config.simulation else None
    assert drive is not None
    assert drive.effective_rate_s is None
    assert drive.intensity_w_m2 == 100.0
    assert drive.dephasing_rad_s == 1.0e7


def test_cascade_state_space_expands_beyond_the_explicit_ladder() -> None:
    root = Path(__file__).resolve().parents[1]
    config = load_transition_config(root / "configs/transitions/three-photon/14P32_13D52.yaml")
    states = states_for_simulation(config)
    assert set(config.ladder).issubset(states)
    assert len(states) > 100
