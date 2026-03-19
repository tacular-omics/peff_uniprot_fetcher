"""Build pefftacular SequenceEntry and FileHeader objects."""

from __future__ import annotations

from collections import Counter

from pefftacular import DatabaseHeader, FileHeader, SequenceEntry

from peff_uniprot_fetcher._fasta import UniProtFastaEntry


def build_entry(fasta_entry: UniProtFastaEntry, annotations: dict) -> SequenceEntry:
    """Build a :class:`SequenceEntry` from a parsed FASTA entry and annotations.

    Parameters
    ----------
    fasta_entry:
        Parsed UniProt FASTA entry.
    annotations:
        Dict returned by :func:`_annotations.features_to_annotations`.
    """
    return SequenceEntry(
        prefix=fasta_entry.db,
        db_unique_id=fasta_entry.accession,
        sequence=fasta_entry.sequence,
        pname=fasta_entry.protein_name or None,
        gname=fasta_entry.gene_name,
        ncbi_tax_id=fasta_entry.tax_id,
        tax_name=fasta_entry.organism,
        length=len(fasta_entry.sequence),
        sv=fasta_entry.sv,
        pe=fasta_entry.pe,
        variant_simple=annotations.get("variant_simple", ()),
        variant_complex=annotations.get("variant_complex", ()),
        mod_res_psi=annotations.get("mod_res_psi", ()),
        mod_res=annotations.get("mod_res", ()),
        processed=annotations.get("processed", ()),
    )


def build_header(entries: list[SequenceEntry], db_version: str | None = None) -> FileHeader:
    """Build a :class:`FileHeader` for the given entries.

    Creates one :class:`DatabaseHeader` per unique prefix (``sp`` /
    ``tr``) found in *entries*.
    """
    counts: Counter[str] = Counter(e.prefix for e in entries)
    databases: list[DatabaseHeader] = []
    for prefix in sorted(counts):
        db_name = "UniProtKB/Swiss-Prot" if prefix == "sp" else "UniProtKB/TrEMBL"
        databases.append(
            DatabaseHeader(
                prefix=prefix,
                db_name=db_name,
                db_version=db_version,
                db_sources=("https://www.uniprot.org",),
                number_of_entries=counts[prefix],
                sequence_type="AA",
            )
        )

    return FileHeader(
        peff_version="1.0",
        databases=tuple(databases),
    )
