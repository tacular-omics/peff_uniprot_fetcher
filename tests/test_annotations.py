"""Tests for GFF feature to PEFF annotation conversion."""

from pefftacular import ModRes, ModResPsi, ModResUnimod, Processed, VariantComplex, VariantSimple
from uniprotptmpy import CrossReference, PtmEntry

from peff_uniprot_fetcher._annotations import features_to_annotations


def _make_ptm(name, psi_mod=None, unimod=None, formula=None, feature_type="MOD_RES", ptm_id=""):  # noqa: E731
    xrefs = []
    if psi_mod:
        xrefs.append(CrossReference("PSI-MOD", psi_mod))
    if unimod is not None:
        xrefs.append(CrossReference("Unimod", str(unimod)))
    return PtmEntry(
        id=ptm_id, name=name, feature_type=feature_type, target="",
        amino_acid_position=None, polypeptide_position=None,
        correction_formula=formula, monoisotopic_mass=None, average_mass=None,
        cellular_location=None, taxonomic_ranges=(), keywords=(),
        cross_references=tuple(xrefs),
    )


PTM_MAP = {
    "Phosphoserine": _make_ptm("Phosphoserine", psi_mod="MOD:00046", unimod=21),
    "Phosphothreonine": _make_ptm("Phosphothreonine", psi_mod="MOD:00047"),
    "UnimodOnly": _make_ptm("UnimodOnly", unimod=340),
    "CustomWithFormula": _make_ptm("CustomWithFormula", formula="C1 H2 O2 S1"),
    "S-palmitoyl cysteine": _make_ptm(
        "S-palmitoyl cysteine", psi_mod="MOD:00111", feature_type="LIPID", ptm_id="PTM-0206",
    ),
    "N-linked (GlcNAc...)": _make_ptm(
        "N-linked (GlcNAc...)", feature_type="CARBOHYD", ptm_id="PTM-0295",
    ),
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
    assert m.name == "M:O-phospho-L-serine"


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
    assert result["mod_res_psi"][0].name == "M:O-phospho-L-serine"


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


def test_modified_residue_unimod_only():
    features = [
        {
            "feature": "Modified residue",
            "start": 400,
            "end": 400,
            "attributes": {"Note": "UnimodOnly"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res_unimod"]) == 1
    m = result["mod_res_unimod"][0]
    assert isinstance(m, ModResUnimod)
    assert m.positions == (400,)
    assert m.accession == "UNIMOD:340"
    assert m.name == "U:Bromo"
    assert result["mod_res_psi"] == ()
    assert result["mod_res"] == ()


def test_modified_residue_custom_with_formula():
    features = [
        {
            "feature": "Modified residue",
            "start": 500,
            "end": 500,
            "attributes": {"Note": "CustomWithFormula"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert result["mod_res_psi"] == ()
    assert result["mod_res_unimod"] == ()
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert isinstance(m, ModRes)
    assert m.accession == "Formula:CH2O2S"
    assert m.name == "CustomWithFormula"


def test_empty_features():
    result = features_to_annotations([], PTM_MAP)
    assert result["variant_simple"] == ()
    assert result["variant_complex"] == ()
    assert result["mod_res_unimod"] == ()
    assert result["mod_res_psi"] == ()
    assert result["mod_res"] == ()
    assert result["processed"] == ()


# -- Glycosylation / Lipidation resolution ---------------------------------


def test_glycosylation_with_ptm_match():
    """Glycosylation whose raw Note matches a PTM map key gets full resolution."""
    features = [
        {
            "feature": "Glycosylation",
            "start": 80,
            "end": 80,
            "attributes": {"Note": "N-linked (GlcNAc...)"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert m.positions == (80,)
    assert m.accession == "PTM-0295"
    assert m.name == "N-linked (GlcNAc...)"


def test_glycosylation_no_ptm_match():
    """Glycosylation with no PTM map match falls back to generic ModRes."""
    features = [
        {
            "feature": "Glycosylation",
            "start": 90,
            "end": 90,
            "attributes": {"Note": "O-linked (Xyl...)"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert m.positions == (90,)
    assert m.accession == ""
    assert m.name == "O-linked"


def test_lipidation_with_ptm_match():
    """Lipidation whose cleaned Note matches a PTM map key gets PSI-MOD resolution."""
    features = [
        {
            "feature": "Lipidation",
            "start": 3,
            "end": 3,
            "attributes": {"Note": "S-palmitoyl cysteine"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    # Should resolve via PSI-MOD cross-reference
    assert len(result["mod_res_psi"]) == 1
    assert result["mod_res_psi"][0].accession == "MOD:00111"
    # Should also get generic ModRes with PTM ID
    assert len(result["mod_res"]) == 1
    assert result["mod_res"][0].accession == "PTM-0206"


def test_lipidation_no_ptm_match():
    """Lipidation with no PTM map match falls back to generic ModRes."""
    features = [
        {
            "feature": "Lipidation",
            "start": 2,
            "end": 2,
            "attributes": {"Note": "GPI-anchor amidated alanine"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert len(result["mod_res"]) == 1
    m = result["mod_res"][0]
    assert m.accession == ""
    assert m.name == "GPI-anchor amidated alanine"


# -- Continue bug fix tests -------------------------------------------------


def test_mod_res_branches_independent():
    """PSI-MOD and UniMod branches don't short-circuit each other."""
    features = [
        {
            "feature": "Modified residue",
            "start": 100,
            "end": 100,
            "attributes": {"Note": "Phosphoserine"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    # Phosphoserine has both PSI-MOD and UniMod xrefs; both should resolve.
    assert len(result["mod_res_psi"]) == 1
    assert len(result["mod_res_unimod"]) == 1
