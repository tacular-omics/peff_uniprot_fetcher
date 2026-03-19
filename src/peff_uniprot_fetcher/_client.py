"""UniProt REST API client."""

import logging

import httpx

log = logging.getLogger(__name__)

UNIPROT_REST_BASE = "https://rest.uniprot.org/uniprotkb"
PTMLIST_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt"


def fetch_entry(accession: str, fmt: str, timeout: float = 30.0) -> str:
    """Fetch a single UniProt entry by accession.

    Parameters
    ----------
    accession:
        UniProt accession (e.g. ``"P12345"``).
    fmt:
        Response format — ``"fasta"`` or ``"gff"``.
    timeout:
        HTTP timeout in seconds.

    Returns
    -------
    str
        The raw response text.
    """
    url = f"{UNIPROT_REST_BASE}/{accession}.{fmt}"
    log.debug("GET %s", url)
    response = httpx.get(url, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    return response.text


def fetch_entries(accessions: list[str], fmt: str, timeout: float = 30.0) -> str:
    """Fetch multiple UniProt entries by accession list.

    Builds a query like ``accession:P12345 OR accession:Q99999`` and hits
    the ``/uniprotkb/search`` endpoint.

    Parameters
    ----------
    accessions:
        List of UniProt accessions.
    fmt:
        Response format — ``"fasta"`` or ``"gff"``.
    timeout:
        HTTP timeout in seconds.

    Returns
    -------
    str
        The concatenated response text.
    """
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    url = f"{UNIPROT_REST_BASE}/search"
    params = {"query": query, "format": fmt, "size": "500"}
    log.info("Fetching %d accession(s) from UniProt (%s)...", len(accessions), fmt)
    response = httpx.get(url, params=params, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    log.debug("Received %d bytes", len(response.content))
    return response.text


def stream_search(query: str, fmt: str, timeout: float = 60.0) -> str:
    """Fetch all search results from UniProt using the stream endpoint.

    Parameters
    ----------
    query:
        UniProt query string (e.g. ``"organism_id:9606 AND reviewed:true"``).
    fmt:
        Response format — ``"fasta"`` or ``"gff"``.
    timeout:
        HTTP timeout in seconds.

    Returns
    -------
    str
        The full response text.
    """
    url = f"{UNIPROT_REST_BASE}/stream"
    params = {"query": query, "format": fmt}
    log.info("Streaming UniProt results for query %r (%s)...", query, fmt)
    response = httpx.get(url, params=params, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    log.info("Received %d bytes", len(response.content))
    return response.text


def fetch_ptmlist(timeout: float = 60.0) -> str:
    """Fetch the UniProt PTM/modification list (ptmlist.txt).

    Returns
    -------
    str
        The raw ptmlist.txt content.
    """
    log.info("Fetching UniProt PTM list...")
    response = httpx.get(PTMLIST_URL, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    log.info("PTM list fetched (%d bytes)", len(response.content))
    return response.text
