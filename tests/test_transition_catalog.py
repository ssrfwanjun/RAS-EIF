from pathlib import Path

from ras_eif.config import load_transition_config


def test_all_legacy_channels_have_explicit_transition_configs() -> None:
    root = Path(__file__).resolve().parents[1]
    configs = sorted((root / "configs/transitions").rglob("*.yaml"))
    assert len(configs) == 31
    loaded = [load_transition_config(path) for path in configs]
    assert {config.identifier for config in loaded} == {path.stem for path in configs}
    assert all(config.expected_thz_ghz > 0.0 for config in loaded)
