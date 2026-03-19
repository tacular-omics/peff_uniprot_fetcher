"""Tests for UniProt FASTA parsing."""

from peff_uniprot_fetcher._fasta import parse_fasta

SAMPLE_FASTA = """\
>sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2
MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYKRDTNSKKM
NLGVGAYRDDNGKPYVLPSVRKAEAQIAAKGLDKEYLPIGGLAEFCRASAELALGENSEV
VKS
>tr|Q99999|Q99999_HUMAN Some protein OS=Homo sapiens OX=9606 PE=4 SV=1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY
"""


def test_parse_fasta_count():
    entries = parse_fasta(SAMPLE_FASTA)
    assert len(entries) == 2


def test_parse_fasta_first_entry():
    entries = parse_fasta(SAMPLE_FASTA)
    e = entries[0]
    assert e.db == "sp"
    assert e.accession == "P12345"
    assert e.entry_name == "AATM_RABIT"
    assert e.protein_name == "Aspartate aminotransferase, mitochondrial"
    assert e.organism == "Oryctolagus cuniculus"
    assert e.tax_id == 9986
    assert e.gene_name == "GOT2"
    assert e.pe == 1
    assert e.sv == 2
    assert e.sequence.startswith("MALLHSARVL")
    assert e.sequence.endswith("VKS")
    assert len(e.sequence) == 123


def test_parse_fasta_trembl_entry():
    entries = parse_fasta(SAMPLE_FASTA)
    e = entries[1]
    assert e.db == "tr"
    assert e.accession == "Q99999"
    assert e.gene_name is None
    assert e.tax_id == 9606
    assert e.pe == 4
    assert e.sv == 1


def test_parse_fasta_empty():
    assert parse_fasta("") == []
    assert parse_fasta("\n\n") == []
