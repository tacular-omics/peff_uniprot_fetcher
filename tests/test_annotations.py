"""Tests for GFF feature to PEFF annotation conversion."""

from pefftacular import ModRes, ModResPsi, Processed, VariantComplex, VariantSimple

from peff_uniprot_fetcher._annotations import features_to_annotations

PTM_MAP = {
    "Phosphoserine": "MOD:00046",
    "Phosphothreonine": "MOD:00047",
}


def test_simple_variant():
    features = [
        {
            "feature": "Natural variant",
            "start": 100,
            "end": 100,
            "attributes": {"Note": "R -> H (dbSNP:rs12345)"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["variant_simple"]) == 1
    v = result["variant_simple"][0]
    assert isinstance(v, VariantSimple)
    assert v.position == 100
    assert v.new_amino_acid == "H"
    assert v.tag == "rs12345"


def test_complex_variant_multichar():
    features = [
        {
            "feature": "Alternative sequence",
            "start": 10,
            "end": 15,
            "attributes": {"Note": "GRSLVK -> AAAAAA"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["variant_complex"]) == 1
    v = result["variant_complex"][0]
    assert isinstance(v, VariantComplex)
    assert v.start_pos == 10
    assert v.end_pos == 15
    assert v.new_sequence == "AAAAAA"


def test_missing_variant():
    features = [
        {
            "feature": "Natural variant",
            "start": 50,
            "end": 55,
            "attributes": {"Note": "Missing"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["variant_complex"]) == 1
    v = result["variant_complex"][0]
    assert v.new_sequence == ""


def test_modified_residue_with_psi():
    features = [
        {
            "feature": "Modified residue",
            "start": 200,
            "end": 200,
            "attributes": {"Note": "Phosphoserine"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res_psi"]) == 1
    m = result["mod_res_psi"][0]
    assert isinstance(m, ModResPsi)
    assert m.positions == (200,)
    assert m.accession == "MOD:00046"
    assert m.name == "Phosphoserine"


def test_modified_residue_without_psi():
    features = [
        {
            "feature": "Modified residue",
            "start": 300,
            "end": 300,
            "attributes": {"Note": "SomeUnknownMod"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert isinstance(m, ModRes)
    assert m.accession == ""
    assert m.name == "SomeUnknownMod"


def test_modified_residue_strips_qualifiers():
    features = [
        {
            "feature": "Modified residue",
            "start": 100,
            "end": 100,
            "attributes": {"Note": "Phosphoserine; alternate"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res_psi"]) == 1
    assert result["mod_res_psi"][0].name == "Phosphoserine"


def test_glycosylation():
    features = [
        {
            "feature": "Glycosylation",
            "start": 50,
            "end": 50,
            "attributes": {"Note": "N-linked (GlcNAc...)"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert m.positions == (50,)


def test_cross_link_spanning():
    features = [
        {
            "feature": "Cross-link",
            "start": 100,
            "end": 200,
            "attributes": {"Note": "Interchain"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert m.positions == (100, 200)


def test_signal_peptide():
    features = [
        {
            "feature": "Signal peptide",
            "start": 1,
            "end": 25,
            "attributes": {},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["processed"]) == 1
    p = result["processed"][0]
    assert isinstance(p, Processed)
    assert p.start_pos == 1
    assert p.end_pos == 25
    assert p.accession == "PEFF:0001001"
    assert p.name == "signal peptide"


def test_chain():
    features = [
        {
            "feature": "Chain",
            "start": 26,
            "end": 430,
            "attributes": {},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["processed"]) == 1
    p = result["processed"][0]
    assert p.accession == "PEFF:0001004"
    assert p.name == "mature protein"


def test_sorting():
    features = [
        {"feature": "Natural variant", "start": 300, "end": 300, "attributes": {"Note": "A -> G"}},
        {"feature": "Natural variant", "start": 100, "end": 100, "attributes": {"Note": "R -> H"}},
        {"feature": "Natural variant", "start": 200, "end": 200, "attributes": {"Note": "K -> E"}},
    ]
    result = features_to_annotations(features, PTM_MAP)
    positions = [v.position for v in result["variant_simple"]]
    assert positions == [100, 200, 300]


def test_empty_features():
    result = features_to_annotations([], PTM_MAP)
    assert result["variant_simple"] == ()
    assert result["variant_complex"] == ()
    assert result["mod_res_psi"] == ()
    assert result["mod_res"] == ()
    assert result["processed"] == ()
