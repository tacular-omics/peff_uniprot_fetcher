"""Tests for PTM map loading and ontology enrichment."""

from uniprotptmpy import CrossReference, PtmEntry

from peff_uniprot_fetcher._ptm import (
    _enrich_from_ontologies,
    get_ptm_map,
    psi_mod_accession,
    unimod_accession,
)


def _make_ptm(**kwargs) -> PtmEntry:
    """Build a minimal PtmEntry with sensible defaults for fields under test."""
    defaults = dict(
        id="PTM-0000",
        name="TestMod",
        feature_type="MOD_RES",
        target="Serine",
        amino_acid_position=None,
        polypeptide_position=None,
        correction_formula=None,
        monoisotopic_mass=None,
        average_mass=None,
        cellular_location=None,
        taxonomic_ranges=(),
        keywords=(),
        cross_references=(),
    )
    defaults.update(kwargs)
    return PtmEntry(**defaults)


# ---------------------------------------------------------------------------
# psi_mod_accession / unimod_accession helpers
# ---------------------------------------------------------------------------


def test_psi_mod_accession_found():
    ptm = _make_ptm(cross_references=(CrossReference("PSI-MOD", "MOD:00046"),))
    assert psi_mod_accession(ptm) == "MOD:00046"


def test_psi_mod_accession_not_found():
    ptm = _make_ptm(cross_references=())
    assert psi_mod_accession(ptm) is None


def test_unimod_accession_found():
    ptm = _make_ptm(cross_references=(CrossReference("Unimod", "21"),))
    assert unimod_accession(ptm) == 21


def test_unimod_accession_not_found():
    ptm = _make_ptm(cross_references=())
    assert unimod_accession(ptm) is None


def test_unimod_accession_multiple_xrefs():
    ptm = _make_ptm(
        cross_references=(
            CrossReference("PSI-MOD", "MOD:00046"),
            CrossReference("Unimod", "21"),
        )
    )
    assert psi_mod_accession(ptm) == "MOD:00046"
    assert unimod_accession(ptm) == 21


# ---------------------------------------------------------------------------
# _enrich_from_ontologies tests
# ---------------------------------------------------------------------------


def test_enrich_unimod_fills_masses():
    """Entry with Unimod xref and no formula should get masses from unimodpy."""
    ptm = _make_ptm(name="Phospho", cross_references=(CrossReference("Unimod", "21"),))
    result = _enrich_from_ontologies({"Phospho": ptm})
    enriched = result["Phospho"]
    assert enriched.monoisotopic_mass is not None
    assert enriched.average_mass is not None
    assert enriched.correction_formula is not None


def test_enrich_already_has_formula_unchanged():
    """Entry that already has a formula must be returned as-is."""
    ptm = _make_ptm(
        name="Phosphoserine",
        correction_formula="H1 O3 P1",
        cross_references=(CrossReference("Unimod", "21"),),
    )
    result = _enrich_from_ontologies({"Phosphoserine": ptm})
    assert result["Phosphoserine"] is ptm


def test_enrich_no_cross_references_unchanged():
    """Entry with no PSI-MOD or Unimod must be returned as-is."""
    ptm = _make_ptm(name="Mystery")
    result = _enrich_from_ontologies({"Mystery": ptm})
    assert result["Mystery"] is ptm


def test_enrich_unrecognised_unimod_unchanged():
    """Entry whose Unimod accession resolves to None is unchanged."""
    ptm = _make_ptm(name="FakeMod", cross_references=(CrossReference("Unimod", "999999"),))
    result = _enrich_from_ontologies({"FakeMod": ptm})
    enriched = result["FakeMod"]
    assert enriched.monoisotopic_mass is None
    assert enriched.average_mass is None
    assert enriched.correction_formula is None


# ---------------------------------------------------------------------------
# get_ptm_map integration test
# ---------------------------------------------------------------------------


def test_get_ptm_map_returns_entries():
    """get_ptm_map should return a non-empty dict of PtmEntry."""
    ptm_map = get_ptm_map()
    assert len(ptm_map) > 0
    first = next(iter(ptm_map.values()))
    assert isinstance(first, PtmEntry)


def test_get_ptm_map_has_phosphoserine():
    """Phosphoserine should be present and have mass data."""
    ptm_map = get_ptm_map()
    ptm = ptm_map["Phosphoserine"]
    assert ptm.monoisotopic_mass is not None
    assert psi_mod_accession(ptm) == "MOD:00046"
