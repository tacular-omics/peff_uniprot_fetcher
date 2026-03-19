"""Tests for PTM list parsing."""

from peff_uniprot_fetcher._ptm import parse_ptmlist

SAMPLE_PTMLIST = """\
------------------------------------------------------------------------
  UniProt Knowledgebase:
  Post-translational modifications
------------------------------------------------------------------------
ID   Phosphoserine
AC   PTM-0001
FT   MOD_RES
DR   PSI-MOD; MOD:00046; O-phospho-L-serine.
//
ID   Phosphothreonine
AC   PTM-0002
FT   MOD_RES
DR   PSI-MOD; MOD:00047; O-phospho-L-threonine.
//
ID   SomeModWithoutPSI
AC   PTM-9999
FT   MOD_RES
//
ID   4-hydroxyproline
AC   PTM-0003
FT   MOD_RES
DR   PSI-MOD; MOD:00039; 4-hydroxy-L-proline.
//
"""


def test_parse_ptmlist_count():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    # Only 3 entries have PSI-MOD; one does not.
    assert len(mapping) == 3


def test_parse_ptmlist_phosphoserine():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert mapping["Phosphoserine"] == "MOD:00046"


def test_parse_ptmlist_phosphothreonine():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert mapping["Phosphothreonine"] == "MOD:00047"


def test_parse_ptmlist_hydroxyproline():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert mapping["4-hydroxyproline"] == "MOD:00039"


def test_parse_ptmlist_missing_psimod():
    mapping = parse_ptmlist(SAMPLE_PTMLIST)
    assert "SomeModWithoutPSI" not in mapping


def test_parse_ptmlist_empty():
    assert parse_ptmlist("") == {}
