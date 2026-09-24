"""Tests for GFF feature to PEFF annotation conversion."""

from pefftacular import ModRes, ModResPsi, ModResUnimod, Processed, VariantComplex, VariantSimple, write_peff
from uniprotptmpy import CrossReference, PtmEntry

from peff_uniprot_fetcher._annotations import features_to_annotations
from peff_uniprot_fetcher._builder import build_entry, build_header
from peff_uniprot_fetcher._fasta import UniProtFastaEntry


def _make_ptm(name, psi_mod=None, unimod=None, formula=None, feature_type="MOD_RES", ptm_id=""):  # noqa: E731
    xrefs = []
    if psi_mod:
        xrefs.append(CrossReference("PSI-MOD", psi_mod))
    if unimod is not None:
        xrefs.append(CrossReference("Unimod", str(unimod)))
    return PtmEntry(
        id=ptm_id,
        name=name,
        feature_type=feature_type,
        target="",
        amino_acid_position=None,
        polypeptide_position=None,
        correction_formula=formula,
        monoisotopic_mass=None,
        average_mass=None,
        cellular_location=None,
        taxonomic_ranges=(),
        keywords=(),
        cross_references=tuple(xrefs),
    )


PTM_MAP = {
    "Phosphoserine": _make_ptm("Phosphoserine", psi_mod="MOD:00046", unimod=21),
    "Phosphothreonine": _make_ptm("Phosphothreonine", psi_mod="MOD:00047"),
    "UnimodOnly": _make_ptm("UnimodOnly", unimod=340),
    "CustomWithFormula": _make_ptm("CustomWithFormula", formula="C1 H2 O2 S1"),
    "S-palmitoyl cysteine": _make_ptm(
        "S-palmitoyl cysteine",
        psi_mod="MOD:00111",
        feature_type="LIPID",
        ptm_id="PTM-0206",
    ),
    "N-linked (GlcNAc...)": _make_ptm(
        "N-linked (GlcNAc...)",
        feature_type="CARBOHYD",
        ptm_id="PTM-0295",
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
    assert m.name == "O-phospho-L-serine"


def test_modified_residue_without_psi():
    """Unknown modifications that don't match ptm_map produce no annotations."""
    features = [
        {
            "feature": "Modified residue",
            "start": 300,
            "end": 300,
            "attributes": {"Note": "SomeUnknownMod"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert result["mod_res"] == ()
    assert result["mod_res_psi"] == ()
    assert result["mod_res_unimod"] == ()


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
    assert result["mod_res_psi"][0].name == "O-phospho-L-serine"


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
    assert m.name == "Bromo"
    assert result["mod_res_psi"] == ()
    # A UNIMOD entry exists, so no generic ModRes (PEFF 1.0 section 3.3.12).
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
    assert m.accession == ""
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
    """Glycosylation with no PTM map match produces no annotations."""
    features = [
        {
            "feature": "Glycosylation",
            "start": 90,
            "end": 90,
            "attributes": {"Note": "O-linked (Xyl...)"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert result["mod_res"] == ()
    assert result["mod_res_psi"] == ()
    assert result["mod_res_unimod"] == ()


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
    # A PSI-MOD entry exists, so no generic ModRes (PEFF 1.0 section 3.3.12).
    assert result["mod_res"] == ()


def test_lipidation_no_ptm_match():
    """Lipidation with no PTM map match produces no annotations."""
    features = [
        {
            "feature": "Lipidation",
            "start": 2,
            "end": 2,
            "attributes": {"Note": "GPI-anchor amidated alanine"},
        }
    ]
    result = features_to_annotations(features, PTM_MAP)
    assert result["mod_res"] == ()
    assert result["mod_res_psi"] == ()
    assert result["mod_res_unimod"] == ()


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
    assert result["mod_res"] == ()


def test_only_known_mass_uses_psimod_delta_mass():
    # MOD:00394 (monoacetylated residue) has a delta mass but no full residue mass;
    # MOD:01928 (N-(L-isoaspartyl)-glycine) has a full residue mass but no delta mass. only_known_mass
    # must judge PSI-MOD by its delta mass, like UNIMOD.
    ptm_map = {
        "DeltaOnly": _make_ptm("DeltaOnly", psi_mod="MOD:00394"),
        "FullOnly": _make_ptm("FullOnly", psi_mod="MOD:01928"),
    }
    features = [
        {"feature": "Modified residue", "start": 1, "end": 1, "attributes": {"Note": "DeltaOnly"}},
        {"feature": "Modified residue", "start": 2, "end": 2, "attributes": {"Note": "FullOnly"}},
    ]
    result = features_to_annotations(features, ptm_map, only_known_mass=True)
    assert [m.accession for m in result["mod_res_psi"]] == ["MOD:00394"]
    result = features_to_annotations(features, ptm_map, only_known_mass=False)
    assert [m.accession for m in result["mod_res_psi"]] == ["MOD:00394", "MOD:01928"]


def test_mod_res_only_without_cv_entry():
    """ModRes is written only when the ptmlist entry maps to neither PSI-MOD nor UNIMOD (section 3.3.12)."""
    ptm_map = {
        "Phosphoserine": _make_ptm("Phosphoserine", psi_mod="MOD:00046", unimod=21, ptm_id="PTM-0253"),
        "NoCv": _make_ptm("NoCv", ptm_id="PTM-9999"),
        "UnknownPsi": _make_ptm("UnknownPsi", psi_mod="MOD:99999", ptm_id="PTM-9998"),
    }
    features = [
        {"feature": "Modified residue", "start": 1, "end": 1, "attributes": {"Note": "Phosphoserine"}},
        {"feature": "Modified residue", "start": 2, "end": 2, "attributes": {"Note": "NoCv"}},
        {"feature": "Modified residue", "start": 3, "end": 3, "attributes": {"Note": "UnknownPsi"}},
    ]
    result = features_to_annotations(features, ptm_map)
    assert [(m.positions, m.accession, m.name) for m in result["mod_res"]] == [
        ((2,), "PTM-9999", "NoCv"),
        ((3,), "PTM-9998", "UnknownPsi"),
    ]
    # only_known_mass dropping a CV entry must not bring back a ModRes for it.
    result = features_to_annotations(features[:1], ptm_map, only_known_mass=True)
    assert result["mod_res"] == ()


def test_peff_output_matches_spec_examples(tmp_path):
    """Written headers match the PEFF 1.0 examples: (100|MOD:00046|O-phospho-L-serine), (100|UNIMOD:21|Phospho)."""
    ptm_map = {"Phosphoserine": _make_ptm("Phosphoserine", psi_mod="MOD:00046", unimod=21, ptm_id="PTM-0253")}
    features = [{"feature": "Modified residue", "start": 100, "end": 100, "attributes": {"Note": "Phosphoserine"}}]
    fasta = UniProtFastaEntry(
        db="sp",
        accession="P12345",
        entry_name="TEST_HUMAN",
        protein_name="Test",
        organism=None,
        tax_id=None,
        gene_name=None,
        pe=None,
        sv=None,
        sequence="S" * 120,
    )
    entry = build_entry(fasta, features_to_annotations(features, ptm_map))
    out = tmp_path / "out.peff"
    write_peff(build_header([entry]), [entry], out)
    header = next(line for line in out.read_text().splitlines() if line.startswith(">sp:P12345"))
    assert r"\ModResPsi=(100|MOD:00046|O-phospho-L-serine)" in header
    assert r"\ModResUnimod=(100|UNIMOD:21|Phospho)" in header
    assert r"\ModRes=" not in header
    assert "M:" not in header and "U:" not in header
