"""peff_uniprot_fetcher -- Generate PEFF files from the UniProt API."""

from __future__ import annotations

import logging
import re
from pathlib import Path

from pefftacular import FileHeader, SequenceEntry, write_peff

from peff_uniprot_fetcher._annotations import features_to_annotations
from peff_uniprot_fetcher._builder import build_entry, build_header
from peff_uniprot_fetcher._client import fetch_entries, stream_search
from peff_uniprot_fetcher._fasta import UniProtFastaEntry, parse_fasta
from peff_uniprot_fetcher._gff import parse_gff
from peff_uniprot_fetcher._ptm import get_ptm_map

log = logging.getLogger(__name__)

_VARIANT_TYPES = {"Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict"}
_PROCESSED_TYPES = {"Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"}

_GFF_MAX_QUERY_LEN = 1800  # conservative limit below UniProt's ~2000-char query cap
_UNIPROT_ACCESSION_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


def _gff_batches(accessions: list[str]) -> list[list[str]]:
    """Split accessions into batches whose OR-query fits within the URL limit."""
    batches: list[list[str]] = []
    batch: list[str] = []
    length = 0
    for acc in accessions:
        term = f"accession:{acc}"
        addition = len(term) + (4 if batch else 0)  # " OR " between terms
        if batch and length + addition > _GFF_MAX_QUERY_LEN:
            batches.append(batch)
            batch = [acc]
            length = len(term)
        else:
            batch.append(acc)
            length += addition
    if batch:
        batches.append(batch)
    return batches


def _fetch_gff_per_accession(accessions: list[str]) -> dict[str, list[dict]]:
    valid = [a for a in accessions if _UNIPROT_ACCESSION_RE.match(a)]
    skipped = len(accessions) - len(valid)
    if skipped:
        log.warning("Skipping %d non-UniProt accession(s) (e.g. contaminants)", skipped)
    total = len(valid)
    all_features: dict[str, list[dict]] = {}
    fetched = 0
    for batch in _gff_batches(valid):
        fetched += len(batch)
        log.info("Fetching GFF annotations: %d / %d", fetched, total)
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        all_features.update(parse_gff(stream_search(query, fmt="gff")))
    return all_features


def _build_entries(
    fasta_entries: list[UniProtFastaEntry],
    all_features: dict[str, list[dict]],
    include_variants: bool,
    include_modifications: bool,
    include_glycosylation: bool,
    include_lipidation: bool,
    include_crosslinks: bool,
    include_processed: bool,
    only_known_mass: bool = False,
) -> tuple[FileHeader, list[SequenceEntry]]:
    if include_modifications:
        ptm_map = get_ptm_map()
        log.info("PTM map loaded (%d entries)", len(ptm_map))
    else:
        ptm_map = {}

    log.info("Building %d PEFF entries...", len(fasta_entries))
    missing_gff = [e.accession for e in fasta_entries if e.accession not in all_features]
    if missing_gff:
        log.warning(
            "%d/%d entries have no GFF features (first: %r). GFF accession keys may not match FASTA accessions.",
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
            elif ft == "Modified residue" and include_modifications:
                filtered.append(feat)
            elif ft == "Glycosylation" and include_glycosylation:
                filtered.append(feat)
            elif ft == "Lipidation" and include_lipidation:
                filtered.append(feat)
            elif ft == "Cross-link" and include_crosslinks:
                filtered.append(feat)
            elif ft in _PROCESSED_TYPES and include_processed:
                filtered.append(feat)

        annotations = features_to_annotations(
            filtered,
            ptm_map,
            only_known_mass=only_known_mass,
        )
        entries.append(build_entry(fasta_entry, annotations))

    log.info("Built %d entries", len(entries))
    return build_header(entries), entries


def fetch_peff(
    accessions: list[str] | None = None,
    query: str | None = None,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_glycosylation: bool = False,
    include_lipidation: bool = False,
    include_crosslinks: bool = False,
    include_processed: bool = True,
    only_known_mass: bool = False,
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
        assert query is not None  # narrowed by the check above
        log.info("Fetching FASTA...")
        fasta_entries = parse_fasta(stream_search(query, fmt="fasta"))
        log.info("Parsed %d sequences from FASTA", len(fasta_entries))
        log.info("Fetching GFF annotations (this may take a while for large result sets)...")
        all_features = parse_gff(stream_search(query, fmt="gff"))
        log.info("Parsed GFF annotations for %d accessions", len(all_features))

    return _build_entries(
        fasta_entries,
        all_features,
        include_variants,
        include_modifications,
        include_glycosylation,
        include_lipidation,
        include_crosslinks,
        include_processed,
        only_known_mass,
    )


def fasta_to_peff(
    fasta: str | Path,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_glycosylation: bool = False,
    include_lipidation: bool = False,
    include_crosslinks: bool = False,
    include_processed: bool = True,
    only_known_mass: bool = False,
) -> tuple[FileHeader, list[SequenceEntry]]:
    """Read a local UniProt FASTA file and annotate from UniProt, returning ``(header, entries)``.

    Sequences come from the local file; GFF annotations are fetched from UniProt per accession.
    """
    log.info("Reading FASTA from %s...", fasta)
    fasta_entries = parse_fasta(Path(fasta).read_text())
    log.info("Parsed %d sequences from FASTA", len(fasta_entries))

    all_features = _fetch_gff_per_accession([e.accession for e in fasta_entries])

    return _build_entries(
        fasta_entries,
        all_features,
        include_variants,
        include_modifications,
        include_glycosylation,
        include_lipidation,
        include_crosslinks,
        include_processed,
        only_known_mass,
    )


def fetch_peff_to_file(
    output: str | Path,
    accessions: list[str] | None = None,
    query: str | None = None,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_glycosylation: bool = False,
    include_lipidation: bool = False,
    include_crosslinks: bool = False,
    include_processed: bool = True,
    only_known_mass: bool = False,
) -> None:
    """Fetch proteins from UniProt and write directly to a PEFF file."""
    header, entries = fetch_peff(
        accessions=accessions,
        query=query,
        include_variants=include_variants,
        include_modifications=include_modifications,
        include_glycosylation=include_glycosylation,
        include_lipidation=include_lipidation,
        include_crosslinks=include_crosslinks,
        include_processed=include_processed,
        only_known_mass=only_known_mass,
    )
    write_peff(header, entries, output)


def fasta_to_peff_file(
    fasta: str | Path,
    output: str | Path,
    include_variants: bool = True,
    include_modifications: bool = True,
    include_glycosylation: bool = False,
    include_lipidation: bool = False,
    include_crosslinks: bool = False,
    include_processed: bool = True,
    only_known_mass: bool = False,
) -> None:
    """Read a local UniProt FASTA file, annotate from UniProt, and write a PEFF file."""
    header, entries = fasta_to_peff(
        fasta=fasta,
        include_variants=include_variants,
        include_modifications=include_modifications,
        include_glycosylation=include_glycosylation,
        include_lipidation=include_lipidation,
        include_crosslinks=include_crosslinks,
        include_processed=include_processed,
        only_known_mass=only_known_mass,
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
