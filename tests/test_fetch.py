"""Tests for accession batching and GFF fetching (UniProt calls are mocked)."""

import logging

import peff_uniprot_fetcher
from peff_uniprot_fetcher import _client

GFF_HEADER = "##gff-version 3\n"


def _gff_line(acc: str) -> str:
    return f"{acc}\tUniProtKB\tModified residue\t15\t15\t.\t.\t.\tNote=Phosphoserine\n"


def _fasta(acc: str) -> str:
    return f">sp|{acc}|TEST_HUMAN Test protein OS=Homo sapiens OX=9606 GN=T PE=1 SV=1\nMSEQ\n"


def _accessions(n: int) -> list[str]:
    return [f"P{i:05d}" for i in range(n)]


def test_accession_regex_accepts_isoforms():
    assert peff_uniprot_fetcher._UNIPROT_ACCESSION_RE.match("P04637")
    assert peff_uniprot_fetcher._UNIPROT_ACCESSION_RE.match("P04637-2")
    assert peff_uniprot_fetcher._UNIPROT_ACCESSION_RE.match("A0A024RBG1-12")
    assert not peff_uniprot_fetcher._UNIPROT_ACCESSION_RE.match("P04637-")
    assert not peff_uniprot_fetcher._UNIPROT_ACCESSION_RE.match("CONT_ALBU_BOVIN")


def test_isoform_accessions_are_queried(monkeypatch):
    queries: list[str] = []

    def fake_stream(query: str, fmt: str, timeout: float = 60.0) -> str:
        queries.append(query)
        return GFF_HEADER + _gff_line("P04637") + _gff_line("P04637-2")

    monkeypatch.setattr(peff_uniprot_fetcher, "stream_search", fake_stream)
    features = peff_uniprot_fetcher._fetch_gff_per_accession(["P04637", "P04637-2"])
    assert "accession:P04637-2" in queries[0]
    assert set(features) == {"P04637", "P04637-2"}


def test_rejected_and_unannotated_isoforms_warn(monkeypatch, caplog):
    monkeypatch.setattr(peff_uniprot_fetcher, "stream_search", lambda q, fmt, timeout=60.0: GFF_HEADER)
    with caplog.at_level(logging.WARNING, logger="peff_uniprot_fetcher"):
        peff_uniprot_fetcher._fetch_gff_per_accession(["P04637-2", "CONT_ALBU_BOVIN"])
    text = caplog.text
    assert "CONT_ALBU_BOVIN" in text  # rejected accessions are named
    assert "P04637-2" in text  # isoforms without features are named


def test_accession_batches_respect_size_and_length():
    accs = _accessions(1200)
    batches = _client.accession_batches(accs, max_query_len=10**9, max_size=500)
    assert [len(b) for b in batches] == [500, 500, 200]
    batches = _client.accession_batches(accs)
    assert [a for b in batches for a in b] == accs
    for b in batches:
        assert len(b) <= _client.SEARCH_PAGE_SIZE
        assert len(" OR ".join(f"accession:{a}" for a in b)) <= _client.MAX_QUERY_LEN


def test_fetch_entries_pages_large_lists(monkeypatch):
    calls: list[dict[str, str]] = []

    def fake_get(url: str, params: dict[str, str] | None, timeout: float) -> str:
        assert params is not None
        calls.append(params)
        accs = [t.removeprefix("accession:") for t in params["query"].split(" OR ")]
        assert len(accs) <= int(params["size"])
        return "".join(_fasta(a) for a in accs)

    monkeypatch.setattr(_client, "_get_text", fake_get)
    accs = _accessions(1200)
    text = _client.fetch_entries(accs, fmt="fasta")
    assert len(calls) > 1
    assert text.count(">sp|") == 1200
    assert all(len(p["query"]) <= _client.MAX_QUERY_LEN for p in calls)


def test_fetch_peff_accessions_over_500(monkeypatch):
    def fake_get(url: str, params: dict[str, str] | None, timeout: float) -> str:
        assert params is not None
        accs = [t.removeprefix("accession:") for t in params["query"].split(" OR ")]
        if params["format"] == "fasta":
            return "".join(_fasta(a) for a in accs[: int(params.get("size", "500"))])
        return GFF_HEADER

    monkeypatch.setattr(_client, "_get_text", fake_get)
    _, entries = peff_uniprot_fetcher.fetch_peff(accessions=_accessions(600), include_modifications=False)
    assert len(entries) == 600
