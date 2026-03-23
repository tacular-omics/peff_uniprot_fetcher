"""Tests for PTM list parsing."""

from peff_uniprot_fetcher._ptm import UniProtPtm, _enrich_from_tacular, parse_ptmlist

SAMPLE_PTMLIST = """\
------------------------------------------------------------------------
  UniProt Knowledgebase:
  Post-translational modifications
------------------------------------------------------------------------
ID   Phosphoserine
AC   PTM-0001
FT   MOD_RES
TG   Serine.
CF   H1 O3 P1
MM   79.966331
MA   79.98
DR   PSI-MOD; MOD:00046; O-phospho-L-serine.
DR   Unimod; 21.
//
ID   Phosphothreonine
AC   PTM-0002
FT   MOD_RES
TG   Threonine.
CF   H1 O3 P1
DR   PSI-MOD; MOD:00047; O-phospho-L-threonine.
//
ID   SomeModWithoutPSI
AC   PTM-9999
FT   MOD_RES
//
ID   4-hydroxyproline
AC   PTM-0003
FT   MOD_RES
TG   Proline.
CF   O1
DR   PSI-MOD; MOD:00039; 4-hydroxy-L-proline.
//
ID   UnimodOnly
AC   PTM-8888
FT   MOD_RES
TG   Tyrosine.
DR   Unimod; 340.
//
"""


def test_parse_ptmlist_count():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert len(mapping) == 5


def test_parse_ptmlist_phosphoserine():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    ptm = mapping["Phosphoserine"]
    assert isinstance(ptm, UniProtPtm)
    assert ptm.psi_mod == "MOD:00046"
    assert ptm.unimod == 21
    assert ptm.formula == "H1 O3 P1"
    assert ptm.proforma_formula == "Formula:H1O3P1"
    assert ptm.mono_mass == 79.966331
    assert ptm.avg_mass == 79.98
    assert ptm.target == "Serine"
    assert ptm.feature_key == "MOD_RES"
    assert ptm.ac == "PTM-0001"


def test_proforma_formula_no_formula():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert mapping["UnimodOnly"].formula is None
    assert mapping["UnimodOnly"].proforma_formula is None


def test_parse_ptmlist_phosphothreonine():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    ptm = mapping["Phosphothreonine"]
    assert ptm.psi_mod == "MOD:00047"
    assert ptm.unimod is None


def test_parse_ptmlist_hydroxyproline():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert mapping["4-hydroxyproline"].psi_mod == "MOD:00039"


def test_parse_ptmlist_no_accessions():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    ptm = mapping["SomeModWithoutPSI"]
    assert ptm.psi_mod is None
    assert ptm.unimod is None


def test_parse_ptmlist_unimod_only():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    ptm = mapping["UnimodOnly"]
    assert ptm.psi_mod is None
    assert ptm.unimod == 340


def test_parse_ptmlist_empty():
    assert parse_ptmlist("") == {}


# ---------------------------------------------------------------------------
# _enrich_from_tacular tests
# ---------------------------------------------------------------------------

def _make_ptm(**kwargs) -> UniProtPtm:
    """Build a minimal UniProtPtm with sensible defaults for fields under test."""
    defaults = dict(
        id="TestMod",
        ac="PTM-0000",
        feature_key="MOD_RES",
        target="Serine",
        formula=None,
        mono_mass=None,
        avg_mass=None,
        psi_mod=None,
        unimod=None,
    )
    defaults.update(kwargs)
    return UniProtPtm(**defaults)


def test_enrich_unimod_fills_masses():
    """Entry with unimod=21 and no formula should get masses from tacular."""
    ptm = _make_ptm(id="Phospho", unimod=21)
    result = _enrich_from_tacular({"Phospho": ptm})
    enriched = result["Phospho"]
    assert enriched.mono_mass is not None
    assert enriched.avg_mass is not None
    assert enriched.formula is not None


def test_enrich_already_has_formula_unchanged():
    """Entry that already has a formula must be returned as-is."""
    ptm = _make_ptm(id="Phosphoserine", formula="H1 O3 P1", unimod=21)
    result = _enrich_from_tacular({"Phosphoserine": ptm})
    assert result["Phosphoserine"] is ptm


def test_enrich_no_cross_references_unchanged():
    """Entry with no psi_mod or unimod must be returned as-is."""
    ptm = _make_ptm(id="Mystery")
    result = _enrich_from_tacular({"Mystery": ptm})
    assert result["Mystery"] is ptm


def test_enrich_unrecognised_unimod_unchanged():
    """Entry whose unimod accession resolves to None in tacular is unchanged."""
    ptm = _make_ptm(id="FakeMod", unimod=999999)
    result = _enrich_from_tacular({"FakeMod": ptm})
    enriched = result["FakeMod"]
    assert enriched.mono_mass is None
    assert enriched.avg_mass is None
    assert enriched.formula is None
