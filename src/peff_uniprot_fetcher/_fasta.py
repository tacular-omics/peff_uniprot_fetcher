"""UniProt FASTA parser."""

from __future__ import annotations

import re
from dataclasses import dataclass

# Matches KEY=value fields in the UniProt FASTA header.
# Keys are 2-character uppercase identifiers: OS, OX, GN, PE, SV.
_KV_PATTERN = re.compile(r"\s(OS|OX|GN|PE|SV)=")


@dataclass
class UniProtFastaEntry:
    """A single parsed UniProt FASTA entry."""

    db: str  # "sp" or "tr"
    accession: str  # e.g. "P12345"
    entry_name: str  # e.g. "PROT_HUMAN"
    protein_name: str
    organism: str | None
    tax_id: int | None
    gene_name: str | None
    pe: int | None
    sv: int | None
    sequence: str


def _parse_header(header: str) -> dict[str, str | None]:
    """Extract key=value fields from a UniProt FASTA header line.

    Returns a dict with keys ``protein_name``, ``OS``, ``OX``, ``GN``,
    ``PE``, ``SV``.  Missing fields are ``None``.
    """
    fields: dict[str, str | None] = {
        "protein_name": None,
        "OS": None,
        "OX": None,
        "GN": None,
        "PE": None,
        "SV": None,
    }

    # Find the first KEY= pattern to split protein name from kv fields.
    first_match = _KV_PATTERN.search(header)
    if first_match is None:
        fields["protein_name"] = header.strip()
        return fields

    fields["protein_name"] = header[: first_match.start()].strip()

    # Walk through all matches to extract key-value pairs.
    matches = list(_KV_PATTERN.finditer(header))
    for i, m in enumerate(matches):
        key = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(header)
        fields[key] = header[start:end].strip()

    return fields


def parse_fasta(text: str) -> list[UniProtFastaEntry]:
    """Parse UniProt FASTA format text into a list of entries."""
    entries: list[UniProtFastaEntry] = []
    current_header: str | None = None
    sequence_parts: list[str] = []

    def _flush() -> None:
        if current_header is None:
            return
        # Header format: >db|accession|entry_name rest_of_header
        pipe1 = current_header.index("|")
        pipe2 = current_header.index("|", pipe1 + 1)
        db = current_header[1:pipe1]  # skip leading ">"
        accession = current_header[pipe1 + 1 : pipe2]
        rest = current_header[pipe2 + 1 :]
        space_idx = rest.find(" ")
        if space_idx == -1:
            entry_name = rest
            header_rest = ""
        else:
            entry_name = rest[:space_idx]
            header_rest = rest[space_idx + 1 :]

        fields = _parse_header(header_rest)
        entries.append(
            UniProtFastaEntry(
                db=db,
                accession=accession,
                entry_name=entry_name,
                protein_name=fields["protein_name"] or "",
                organism=fields["OS"],
                tax_id=int(fields["OX"]) if fields["OX"] else None,
                gene_name=fields["GN"],
                pe=int(fields["PE"]) if fields["PE"] else None,
                sv=int(fields["SV"]) if fields["SV"] else None,
                sequence="".join(sequence_parts),
            )
        )

    for line in text.splitlines():
        if line.startswith(">"):
            _flush()
            current_header = line
            sequence_parts = []
        else:
            sequence_parts.append(line.strip())

    _flush()
    return entries
