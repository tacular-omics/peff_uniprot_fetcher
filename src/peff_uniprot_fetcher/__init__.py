"""peff_uniprot_fetcher -- Generate PEFF files from the UniProt API."""

from __future__ import annotations

import logging
import re
from pathlib import Path

from pefftacular import FileHeader, SequenceEntry, write_peff

from peff_uniprot_fetcher._annotations import features_to_annotations
from peff_uniprot_fetcher._builder import build_entry, build_header
from peff_uniprot_fetcher._client import accession_batches, fetch_entries, stream_search
from peff_uniprot_fetcher._config import AnnotationConfig
from peff_uniprot_fetcher._fasta import UniProtFastaEntry, parse_fasta
from peff_uniprot_fetcher._gff import parse_gff
from peff_uniprot_fetcher._ptm import get_ptm_map

log = logging.getLogger(__name__)

_VARIANT_TYPES = {"Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict"}
_PROCESSED_TYPES = {"Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"}

# Primary accessions (https://www.uniprot.org/help/accession_numbers), optionally with an
# isoform suffix such as "-2".
_UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$"
)
_MAX_LOGGED = 10  # how many example accessions to name in warnings


def _examples(accessions: list[str]) -> str:
    shown = ", ".join(accessions[:_MAX_LOGGED])
    return shown + (f", ... (+{len(accessions) - _MAX_LOGGED} more)" if len(accessions) > _MAX_LOGGED else "")


def _fetch_gff_per_accession(accessions: list[str]) -> dict[str, list[dict]]:
    valid = [a for a in accessions if _UNIPROT_ACCESSION_RE.match(a)]
    rejected = [a for a in accessions if not _UNIPROT_ACCESSION_RE.match(a)]
    if rejected:
        log.warning(
            "%d accession(s) are not UniProt accessions and get NO annotations (e.g. contaminants): %s",
            len(rejected),
            _examples(rejected),
        )
    total = len(valid)
    all_features: dict[str, list[dict]] = {}
    fetched = 0
    for batch in accession_batches(valid):
        fetched += len(batch)
        log.info("Fetching GFF annotations: %d / %d", fetched, total)
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        all_features.update(parse_gff(stream_search(query, fmt="gff")))
    bare_isoforms = [a for a in valid if "-" in a and a not in all_features]
    if bare_isoforms:
        log.warning(
            "%d isoform accession(s) got NO annotations: UniProt publishes GFF features only for "
            "canonical sequences, and canonical positions do not apply to other isoforms: %s",
            len(bare_isoforms),
            _examples(bare_isoforms),
        )
    return all_features


def _build_entries(
    fasta_entries: list[UniProtFastaEntry],
    all_features: dict[str, list[dict]],
    cfg: AnnotationConfig,
) -> tuple[FileHeader, list[SequenceEntry]]:
    if cfg.include_modifications:
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
            if ft in _VARIANT_TYPES and cfg.include_variants:
                filtered.append(feat)
            elif ft == "Modified residue" and cfg.include_modifications:
                filtered.append(feat)
            elif ft == "Glycosylation" and cfg.include_glycosylation:
                filtered.append(feat)
            elif ft == "Lipidation" and cfg.include_lipidation:
                filtered.append(feat)
            elif ft == "Cross-link" and cfg.include_crosslinks:
                filtered.append(feat)
            elif ft in _PROCESSED_TYPES and cfg.include_processed:
                filtered.append(feat)

        annotations = features_to_annotations(
            filtered,
            ptm_map,
            only_known_mass=cfg.only_known_mass,
        )
        entries.append(build_entry(fasta_entry, annotations))

    log.info("Built %d entries", len(entries))
    return build_header(entries), entries


def fetch_peff(
    accessions: list[str] | None = None,
    query: str | None = None,
    cfg: AnnotationConfig | None = None,
    **kwargs: bool,
) -> tuple[FileHeader, list[SequenceEntry]]:
    """Fetch proteins from UniProt and return ``(header, entries)`` for writing as PEFF.

    Either *accessions* or *query* must be provided (not both).

    Annotation behaviour is controlled via *cfg* (:class:`AnnotationConfig`).
    For convenience, individual flags (e.g. ``include_variants=False``) may be
    passed as keyword arguments instead; they are forwarded to the
    :class:`AnnotationConfig` constructor.
    """
    if cfg is None:
        cfg = AnnotationConfig(**kwargs)
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

    return _build_entries(fasta_entries, all_features, cfg)


def fasta_to_peff(
    fasta: str | Path,
    cfg: AnnotationConfig | None = None,
    **kwargs: bool,
) -> tuple[FileHeader, list[SequenceEntry]]:
    """Read a local UniProt FASTA file and annotate from UniProt, returning ``(header, entries)``.

    Sequences come from the local file; GFF annotations are fetched from UniProt per accession.
    See :func:`fetch_peff` for annotation options.
    """
    if cfg is None:
        cfg = AnnotationConfig(**kwargs)
    log.info("Reading FASTA from %s...", fasta)
    fasta_entries = parse_fasta(Path(fasta).read_text())
    log.info("Parsed %d sequences from FASTA", len(fasta_entries))

    all_features = _fetch_gff_per_accession([e.accession for e in fasta_entries])

    return _build_entries(fasta_entries, all_features, cfg)


def fetch_peff_to_file(
    output: str | Path,
    accessions: list[str] | None = None,
    query: str | None = None,
    cfg: AnnotationConfig | None = None,
    **kwargs: bool,
) -> None:
    """Fetch proteins from UniProt and write directly to a PEFF file."""
    if cfg is None:
        cfg = AnnotationConfig(**kwargs)
    header, entries = fetch_peff(accessions=accessions, query=query, cfg=cfg)
    write_peff(header, entries, output)


def fasta_to_peff_file(
    fasta: str | Path,
    output: str | Path,
    cfg: AnnotationConfig | None = None,
    **kwargs: bool,
) -> None:
    """Read a local UniProt FASTA file, annotate from UniProt, and write a PEFF file."""
    if cfg is None:
        cfg = AnnotationConfig(**kwargs)
    header, entries = fasta_to_peff(fasta=fasta, cfg=cfg)
    write_peff(header, entries, output)


__all__ = [
    "AnnotationConfig",
    "FileHeader",
    "SequenceEntry",
    "fasta_to_peff",
    "fasta_to_peff_file",
    "fetch_peff",
    "fetch_peff_to_file",
    "write_peff",
]
