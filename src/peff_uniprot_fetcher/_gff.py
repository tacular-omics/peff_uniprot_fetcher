"""UniProt GFF3 feature parser."""

from __future__ import annotations

import logging
from urllib.parse import unquote

log = logging.getLogger(__name__)

# Feature types we extract from UniProt GFF.
VARIANT_TYPES = frozenset(
    {
        "Natural variant",
        "Mutagenesis",
        "Alternative sequence",
        "Sequence conflict",
    }
)

MODIFICATION_TYPES = frozenset(
    {
        "Modified residue",
        "Glycosylation",
        "Lipidation",
        "Cross-link",
    }
)

PROCESSED_TYPES = frozenset(
    {
        "Signal peptide",
        "Transit peptide",
        "Propeptide",
        "Chain",
        "Peptide",
    }
)

_ALL_TYPES = VARIANT_TYPES | MODIFICATION_TYPES | PROCESSED_TYPES


def _parse_attributes(raw: str) -> dict[str, str]:
    """Parse a GFF3 attributes column (column 9) into a dict.

    URL-decodes the entire string first, then splits on ``;`` only at
    parenthesis depth 0 so that semicolons inside values like
    ``N-linked (GlcNAc...; detail)`` are preserved correctly.
    """
    decoded = unquote(raw)
    attrs: dict[str, str] = {}
    depth = 0
    part_start = 0

    def _store(segment: str) -> None:
        segment = segment.strip()
        if "=" in segment:
            key, _, value = segment.partition("=")
            attrs[key.strip()] = value.strip()

    for i, ch in enumerate(decoded):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == ";" and depth == 0:
            _store(decoded[part_start:i])
            part_start = i + 1

    _store(decoded[part_start:])
    return attrs


def parse_gff(gff_text: str) -> dict[str, list[dict]]:
    """Parse GFF3 text and return features grouped by accession.

    Parameters
    ----------
    gff_text:
        Raw GFF3 text from UniProt.

    Returns
    -------
    dict[str, list[dict]]
        Mapping of accession to list of feature dicts. Each dict has
        keys ``feature``, ``start``, ``end``, and ``attributes``.
    """
    result: dict[str, list[dict]] = {}
    seen_types: set[str] = set()
    data_lines = 0

    for line in gff_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue

        parts = line.split("\t")
        if len(parts) < 9:
            log.debug("Skipping line with %d columns: %r", len(parts), line[:120])
            continue

        data_lines += 1
        feature_type = parts[2]
        seen_types.add(feature_type)
        if feature_type not in _ALL_TYPES:
            continue

        # UniProt GFF column 1 is always the bare accession (e.g. "P12345"),
        # both in per-entry and stream responses.
        accession = parts[0]

        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue

        attributes = _parse_attributes(parts[8])

        feature_dict = {
            "feature": feature_type,
            "start": start,
            "end": end,
            "attributes": attributes,
        }

        result.setdefault(accession, []).append(feature_dict)

    total = sum(len(v) for v in result.values())
    log.info("GFF parsed: %d features across %d entries (saw %d data lines)", total, len(result), data_lines)
    if seen_types:
        log.info("Feature types in GFF: %s", sorted(seen_types))

    return result
