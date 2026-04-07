"""Browser (Pyodide) entry point for peff_uniprot_fetcher.

This module exposes a single function, :func:`generate_peff_string`, that the
GitHub Pages + Pyodide frontend calls. It wraps :func:`fetch_peff` and writes
the result into an in-memory string buffer, so no filesystem is needed.
"""

from __future__ import annotations

import io

from pefftacular import write_peff

from peff_uniprot_fetcher import fetch_peff
from peff_uniprot_fetcher._config import AnnotationConfig


def generate_peff_string(organism_id: int | str, reviewed: bool = True) -> str:
    """Generate a PEFF file for *organism_id* as a single string.

    Parameters
    ----------
    organism_id:
        NCBI taxonomy ID (e.g. ``9606`` for human, ``83333`` for E. coli).
    reviewed:
        If ``True`` (default), restrict to reviewed (Swiss-Prot) entries.
    """
    reviewed_clause = " AND reviewed:true" if reviewed else ""
    query = f"organism_id:{organism_id}{reviewed_clause}"
    header, entries = fetch_peff(query=query, cfg=AnnotationConfig())
    buf = io.StringIO()
    write_peff(header, entries, buf)
    return buf.getvalue()
