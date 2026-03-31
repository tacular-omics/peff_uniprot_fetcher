"""Tests for CLI argument parsing and AnnotationConfig construction."""

from __future__ import annotations

import argparse

from peff_uniprot_fetcher._cli import _annotation_config, _annotation_flags
from peff_uniprot_fetcher._config import AnnotationConfig


def _parse(args: list[str]) -> argparse.Namespace:
    """Build a minimal parser with annotation flags and parse *args*."""
    parser = argparse.ArgumentParser()
    _annotation_flags(parser)
    return parser.parse_args(args)


def test_default_flags():
    cfg = _annotation_config(_parse([]))
    assert cfg == AnnotationConfig()
    assert cfg.include_variants is True
    assert cfg.include_modifications is True
    assert cfg.include_processed is True
    assert cfg.include_glycosylation is False
    assert cfg.include_lipidation is False
    assert cfg.include_crosslinks is False
    assert cfg.only_known_mass is False


def test_disable_flags():
    cfg = _annotation_config(_parse(["--no-variants", "--no-modifications", "--no-processed"]))
    assert cfg.include_variants is False
    assert cfg.include_modifications is False
    assert cfg.include_processed is False


def test_enable_optional_flags():
    cfg = _annotation_config(_parse(["--glycosylation", "--lipidation", "--crosslinks", "--only-known-mass"]))
    assert cfg.include_glycosylation is True
    assert cfg.include_lipidation is True
    assert cfg.include_crosslinks is True
    assert cfg.only_known_mass is True
