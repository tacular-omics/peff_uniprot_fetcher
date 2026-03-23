"""Tests for PTM list parsing."""

from peff_uniprot_fetcher._ptm import UniProtPtm, parse_ptmlist

SAMPLE_PTMLIST = """\
------------------------------------------------------------------------
  UniProt Knowledgebase:
  Post-translational modifications
------------------------------------------------------------------------
ID   Phosphoserine
AC   PTM-0001
FT   MOD_RES
TG   Serine.
MM   79.966331
MA   79.98
DR   PSI-MOD; MOD:00046; O-phospho-L-serine.
DR   Unimod; 21.
//
ID   Phosphothreonine
AC   PTM-0002
FT   MOD_RES
TG   Threonine.
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
    assert ptm.mono_mass == 79.966331
    assert ptm.avg_mass == 79.98
    assert ptm.target == "Serine"
    assert ptm.feature_key == "MOD_RES"
    assert ptm.ac == "PTM-0001"


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
