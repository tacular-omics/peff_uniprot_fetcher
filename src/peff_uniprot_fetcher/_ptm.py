"""PTM name to PSI-MOD accession mapping."""

from peff_uniprot_fetcher._client import fetch_ptmlist

_PTM_CACHE: dict[str, str] | None = None


def parse_ptmlist(text: str) -> dict[str, str]:
    """Parse UniProt ptmlist.txt text and return ``{ptm_name: psi_mod_accession}``.

    Records are separated by ``//`` lines.  Each record has an ``ID``
    line and optionally a ``DR   PSI-MOD; MOD:XXXXX; ...`` line.
    """
    mapping: dict[str, str] = {}
    current_id: str | None = None
    current_psimod: str | None = None

    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("ID"):
            current_id = stripped[5:].strip()
        elif stripped.startswith("DR   PSI-MOD;"):
            current_psimod = stripped.split(";")[1].strip().rstrip(".")
        elif stripped == "//":
            if current_id and current_psimod:
                mapping[current_id] = current_psimod
            current_id = None
            current_psimod = None

    return mapping


def get_ptm_map() -> dict[str, str]:
    """Fetch and return the PTM map, caching after the first call."""
    global _PTM_CACHE  # noqa: PLW0603
    if _PTM_CACHE is None:
        _PTM_CACHE = parse_ptmlist(fetch_ptmlist())
    return _PTM_CACHE
