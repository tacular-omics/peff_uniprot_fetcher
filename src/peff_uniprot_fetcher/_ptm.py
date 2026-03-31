"""PTM name to PtmEntry mapping with ontology enrichment."""

from __future__ import annotations

from dataclasses import replace

import psimodpy
import unimodpy
import uniprotptmpy
from uniprotptmpy import PtmEntry

_PTM_CACHE: dict[str, PtmEntry] | None = None

_psimod_db: psimodpy.PsiModDatabase | None = None
_unimod_db: unimodpy.UnimodDatabase | None = None


def _get_psimod_db() -> psimodpy.PsiModDatabase:
    global _psimod_db  # noqa: PLW0603
    if _psimod_db is None:
        _psimod_db = psimodpy.load()
    return _psimod_db


def _get_unimod_db() -> unimodpy.UnimodDatabase:
    global _unimod_db  # noqa: PLW0603
    if _unimod_db is None:
        _unimod_db = unimodpy.load()
    return _unimod_db


def get_psimod(id: str) -> psimodpy.PsiModEntry | None:
    """Get PSI-MOD info for a given accession (e.g. ``"MOD:00046"``), or ``None`` if not found."""
    return _get_psimod_db().get_by_id(id)


def get_unimod(id: int) -> unimodpy.UnimodEntry | None:
    """Get UniMod info for a given numeric accession (e.g. 35), or ``None`` if not found."""
    return _get_unimod_db().get_by_id(id)


def psi_mod_accession(entry: PtmEntry) -> str | None:
    """Extract PSI-MOD accession (e.g. ``"MOD:00046"``) from cross-references, or ``None``."""
    for xref in entry.cross_references:
        if xref.database == "PSI-MOD":
            return xref.accession
    return None


def unimod_accession(entry: PtmEntry) -> int | None:
    """Extract numeric UniMod accession from cross-references, or ``None``."""
    for xref in entry.cross_references:
        if xref.database == "Unimod":
            try:
                return int(xref.accession)
            except ValueError:
                pass
    return None


def _enrich_from_ontologies(
    ptm_map: dict[str, PtmEntry],
) -> dict[str, PtmEntry]:
    """Fill correction_formula / monoisotopic_mass / average_mass from unimodpy/psimodpy
    for entries that have a PSI-MOD or UniMod cross-reference but no mass data from UniProt."""
    enriched: dict[str, PtmEntry] = {}
    for name, ptm in ptm_map.items():
        unimod_id = unimod_accession(ptm)
        psi_mod_id = psi_mod_accession(ptm)

        if ptm.correction_formula is not None or (unimod_id is None and psi_mod_id is None):
            enriched[name] = ptm
            continue

        # Prefer UniMod (more standardised masses), fall back to PSI-MOD
        formula: str | None = None
        mono_mass: float | None = None
        avg_mass: float | None = None

        if unimod_id is not None:
            info = _get_unimod_db().get_by_id(unimod_id)
            if info is not None:
                if info.dict_composition:
                    formula = " ".join(f"{el}{cnt}" for el, cnt in sorted(info.dict_composition.items()))
                mono_mass = info.delta_mono_mass
                avg_mass = info.delta_avge_mass

        if formula is None and mono_mass is None and psi_mod_id is not None:
            try:
                psi_num = int(psi_mod_id.split(":")[1])
            except (IndexError, ValueError):
                psi_num = None
            if psi_num is not None:
                info = _get_psimod_db().get_by_id(psi_num)
                if info is not None:
                    if info.dict_diff_formula:
                        formula = " ".join(f"{el}{cnt}" for el, cnt in sorted(info.dict_diff_formula.items()))
                    mono_mass = info.diff_mono
                    avg_mass = info.diff_avg

        if formula is None and mono_mass is None:
            enriched[name] = ptm
            continue

        enriched[name] = replace(
            ptm,
            correction_formula=formula if formula is not None else ptm.correction_formula,
            monoisotopic_mass=mono_mass if ptm.monoisotopic_mass is None else ptm.monoisotopic_mass,
            average_mass=avg_mass if ptm.average_mass is None else ptm.average_mass,
        )
    return enriched


def get_ptm_map() -> dict[str, PtmEntry]:
    """Load and return the PTM map, caching after the first call."""
    global _PTM_CACHE  # noqa: PLW0603
    if _PTM_CACHE is None:
        db = uniprotptmpy.load()
        ptm_map = {entry.name: entry for entry in db}
        _PTM_CACHE = _enrich_from_ontologies(ptm_map)
    return _PTM_CACHE
