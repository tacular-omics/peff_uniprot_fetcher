"""UniProt REST API client.

Uses ``httpx`` when running on a normal Python interpreter. When running under
Pyodide (``sys.platform == "emscripten"``), httpx's socket-based transport is
unavailable, so we route through ``pyodide.http.pyfetch`` instead. This lets
the same package power the CLI and a browser (Pyodide) frontend.
"""

import logging
import sys
from urllib.parse import urlencode

log = logging.getLogger(__name__)

UNIPROT_REST_BASE = "https://rest.uniprot.org/uniprotkb"

_IS_PYODIDE = sys.platform == "emscripten"

if not _IS_PYODIDE:
    import httpx


def _get_text(url: str, params: dict[str, str] | None, timeout: float) -> str:
    """GET *url* with *params* and return the response body as text.

    Dispatches to ``httpx`` natively and to ``pyodide.http.pyfetch`` under
    Pyodide. The Pyodide branch uses a synchronous XHR fallback via
    ``pyfetch`` awaited on the event loop, which works from a Web Worker.
    """
    if _IS_PYODIDE:
        # Use synchronous XHR via JS. Sync XHR is permitted inside Web
        # Workers (which is where we run Pyodide on the site), and it keeps
        # the rest of ``peff_uniprot_fetcher`` synchronous. Do NOT call this
        # from the main thread.
        from js import XMLHttpRequest  # type: ignore[import-not-found]

        full_url = url if not params else f"{url}?{urlencode(params)}"
        xhr = XMLHttpRequest.new()
        xhr.open("GET", full_url, False)  # False = synchronous
        xhr.send(None)
        if xhr.status >= 400:
            raise RuntimeError(f"UniProt HTTP {xhr.status} for {full_url}")
        return xhr.responseText

    response = httpx.get(url, params=params, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    return response.text


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
    return _get_text(url, None, timeout)


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
    text = _get_text(url, params, timeout)
    log.debug("Received %d bytes", len(text))
    return text


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
    text = _get_text(url, params, timeout)
    log.info("Received %d bytes", len(text))
    return text
