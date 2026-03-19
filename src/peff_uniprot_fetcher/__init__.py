"""peff_uniprot_fetcher -- Generate PEFF files from the UniProt API."""

from __future__ import annotations

import logging
from pathlib import Path

from pefftacular import FileHeader, SequenceEntry, write_peff

from peff_uniprot_fetcher._annotations import features_to_annotations
from peff_uniprot_fetcher._builder import build_entry, build_header
from peff_uniprot_fetcher._client import fetch_entries, fetch_entry, stream_search
from peff_uniprot_fetcher._fasta import UniProtFastaEntry, parse_fasta
from peff_uniprot_fetcher._gff import parse_gff
from peff_uniprot_fetcher._ptm import get_ptm_map

log = logging.getLogger(__name__)

_VARIANT_TYPES = {"Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict"}
_MOD_TYPES = {"Modified residue", "Glycosylation", "Lipidation", "Cross-link"}
_PROCESSED_TYPES = {"Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"}

_GFF_LOG_INTERVAL = 100  # log progress every N accessions during per-entry GFF fetching


def _fetch_gff_per_accession(accessions: list[str]) -> dict[str, list[dict]]:
    total = len(accessions)
    all_features: dict[str, list[dict]] = {}
    for i, acc in enumerate(accessions, 1):
        if i == 1 or i % _GFF_LOG_INTERVAL == 0 or i == total:
            log.info("Fetching GFF annotations: %d / %d", i, total)
        all_features.update(parse_gff(fetch_entry(acc, fmt="gff")))
    return all_features


def _build_entries(
    fasta_entries: list[UniProtFastaEntry],
    all_features: dict[str, list[dict]],
    include_variants: bool,
    include_modifications: bool,
    include_processed: bool,
) -> tuple[FileHeader, list[SequenceEntry]]:
    if include_modifications:
        ptm_map: dict[str, str] = get_ptm_map()
        log.info("PTM map loaded (%d entries)", len(ptm_map))
    else:
        ptm_map = {}

    log.info("Building %d PEFF entries...", len(fasta_entries))
    missing_gff = [e.accession for e in fasta_entries if e.accession not in all_features]
    if missing_gff:
        log.warning(
            "%d/%d entries have no GFF features (first: %r). "
            "GFF accession keys may not match FASTA accessions.",
            len(missing_gff),
            len(fasta_entries),
            missing_gff[0],
        )
    entries: list[SequenceEntry] = []
    for fasta_entry in fasta_entries:
        raw_features = all_features.get(fasta_entry.accession, [])

        filtered: list[dict] = []
        for feat in raw_features:
            ft = feat["feature"]
            if ft in _VARIANT_TYPES and include_variants:
                filtered.append(feat)
            elif ft in _MOD_TYPES and include_modifications:
                filtered.append(feat)
            elif ft in _PROCESSED_TYPES and include_processed:
                filtered.append(feat)

        annotations = features_to_annotations(filtered, ptm_map)
        entries.append(build_entry(fasta_entry, annotations))

    log.info("Built %d entries", len(entries))
    return build_header(entries), entries


def fetch_peff(
    accessions: list[str] | None = None,
    query: str | None = None,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_processed: bool = True,
) -> tuple[FileHeader, list[SequenceEntry]]:
    """Fetch proteins from UniProt and return ``(header, entries)`` for writing as PEFF.

    Either *accessions* or *query* must be provided (not both).
    """
    if (accessions is None) == (query is None):
        raise ValueError("Provide exactly one of 'accessions' or 'query'.")

    if accessions is not None:
        fasta_entries = parse_fasta(fetch_entries(accessions, fmt="fasta"))
        log.info("Parsed %d sequences from FASTA", len(fasta_entries))
        all_features = _fetch_gff_per_accession(accessions)
    else:
        log.info("Fetching FASTA...")
        fasta_entries = parse_fasta(stream_search(query, fmt="fasta"))  # type: ignore[arg-type]
        log.info("Parsed %d sequences from FASTA", len(fasta_entries))
        log.info("Fetching GFF annotations (this may take a while for large result sets)...")
        all_features = parse_gff(stream_search(query, fmt="gff"))  # type: ignore[arg-type]
        log.info("Parsed GFF annotations for %d accessions", len(all_features))

    return _build_entries(fasta_entries, all_features, include_variants, include_modifications, include_processed)


def fasta_to_peff(
    fasta: str | Path,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_processed: bool = True,
) -> tuple[FileHeader, list[SequenceEntry]]:
    """Read a local UniProt FASTA file and annotate from UniProt, returning ``(header, entries)``.

    Sequences come from the local file; GFF annotations are fetched from UniProt per accession.
    """
    log.info("Reading FASTA from %s...", fasta)
    fasta_entries = parse_fasta(Path(fasta).read_text())
    log.info("Parsed %d sequences from FASTA", len(fasta_entries))

    all_features = _fetch_gff_per_accession([e.accession for e in fasta_entries])

    return _build_entries(fasta_entries, all_features, include_variants, include_modifications, include_processed)


def fetch_peff_to_file(
    output: str | Path,
    accessions: list[str] | None = None,
    query: str | None = None,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_processed: bool = True,
) -> None:
    """Fetch proteins from UniProt and write directly to a PEFF file."""
    header, entries = fetch_peff(
        accessions=accessions,
        query=query,
        include_variants=include_variants,
        include_modifications=include_modifications,
        include_processed=include_processed,
    )
    write_peff(header, entries, output)


def fasta_to_peff_file(
    fasta: str | Path,
    output: str | Path,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_processed: bool = True,
) -> None:
    """Read a local UniProt FASTA file, annotate from UniProt, and write a PEFF file."""
    header, entries = fasta_to_peff(
        fasta=fasta,
        include_variants=include_variants,
        include_modifications=include_modifications,
        include_processed=include_processed,
    )
    write_peff(header, entries, output)


__all__ = [
    "FileHeader",
    "SequenceEntry",
    "fasta_to_peff",
    "fasta_to_peff_file",
    "fetch_peff",
    "fetch_peff_to_file",
    "write_peff",
]
