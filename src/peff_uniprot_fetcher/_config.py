"""Annotation configuration for PEFF generation."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class AnnotationConfig:
    """Controls which annotation types are included in PEFF output."""

    include_variants: bool = True
    include_modifications: bool = True
    include_glycosylation: bool = False
    include_lipidation: bool = False
    include_crosslinks: bool = False
    include_processed: bool = True
    only_known_mass: bool = False
