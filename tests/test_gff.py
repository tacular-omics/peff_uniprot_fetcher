"""Tests for UniProt GFF3 parsing."""

from peff_uniprot_fetcher._gff import parse_gff

SAMPLE_GFF = """\
##gff-version 3
##sequence-region P12345 1 430
P12345\tUniProtKB\tNatural variant\t100\t100\t.\t.\t.\tID=VAR_012345;Note=R -> H (dbSNP:rs12345)
P12345\tUniProtKB\tModified residue\t200\t200\t.\t.\t.\tNote=Phosphoserine
P12345\tUniProtKB\tSignal peptide\t1\t25\t.\t.\t.\tNote=Signal peptide
P12345\tUniProtKB\tChain\t26\t430\t.\t.\t.\tNote=Mature protein
P12345\tUniProtKB\tDisulfide bond\t150\t250\t.\t.\t.\tNote=Interchain
Q99999\tUniProtKB\tGlycosylation\t50\t50\t.\t.\t.\tNote=N-linked (GlcNAc...)
Q99999\tUniProtKB\tAlternative sequence\t10\t15\t.\t.\t.\tNote=GRSLVK -> AAAAAA; in isoform 2
"""


def test_parse_gff_grouping():
    result = parse_gff(SAMPLE_GFF)
    assert "P12345" in result
    assert "Q99999" in result


def test_parse_gff_p12345_features():
    result = parse_gff(SAMPLE_GFF)
    features = result["P12345"]
    # Natural variant, Modified residue, Signal peptide, Chain (Disulfide bond excluded)
    assert len(features) == 4
    types = [f["feature"] for f in features]
    assert "Natural variant" in types
    assert "Modified residue" in types
    assert "Signal peptide" in types
    assert "Chain" in types
    assert "Disulfide bond" not in types


def test_parse_gff_positions():
    result = parse_gff(SAMPLE_GFF)
    variant = [f for f in result["P12345"] if f["feature"] == "Natural variant"][0]
    assert variant["start"] == 100
    assert variant["end"] == 100


def test_parse_gff_attributes():
    result = parse_gff(SAMPLE_GFF)
    variant = [f for f in result["P12345"] if f["feature"] == "Natural variant"][0]
    assert "Note" in variant["attributes"]
    assert "R -> H" in variant["attributes"]["Note"]


def test_parse_gff_q99999():
    result = parse_gff(SAMPLE_GFF)
    features = result["Q99999"]
    assert len(features) == 2
    types = [f["feature"] for f in features]
    assert "Glycosylation" in types
    assert "Alternative sequence" in types


def test_parse_gff_skips_comments():
    gff = "##gff-version 3\n# a comment\n"
    result = parse_gff(gff)
    assert result == {}


def test_parse_gff_empty():
    assert parse_gff("") == {}
