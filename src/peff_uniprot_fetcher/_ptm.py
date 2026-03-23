"""PTM name to UniProt PTM entry mapping."""

from __future__ import annotations

from dataclasses import dataclass

from peff_uniprot_fetcher._client import fetch_ptmlist

_PTM_CACHE: dict[str, UniProtPtm] | None = None


@dataclass(frozen=True)
class UniProtPtm:
    """A single entry from UniProt's ptmlist.txt.

    Attributes
    ----------
    id:
        PTM identifier as it appears in GFF ``Note`` fields (e.g. ``"Phosphoserine"``).
    ac:
        UniProt internal accession (e.g. ``"PTM-0253"``).
    feature_key:
        UniProt feature type: ``MOD_RES``, ``CROSSLNK``, ``LIPID``, ``CARBOHYD``, etc.
    target:
        Amino acid target(s) (e.g. ``"Serine"``).
    mono_mass:
        Monoisotopic mass difference in Da, or ``None`` if not listed.
    avg_mass:
        Average mass difference in Da, or ``None`` if not listed.
    psi_mod:
        PSI-MOD accession (e.g. ``"MOD:00046"``), or ``None`` if not cross-referenced.
    unimod:
        UniMod accession number (e.g. ``21``), or ``None`` if not cross-referenced.
    """

    id: str
    ac: str
    feature_key: str
    target: str
    mono_mass: float | None
    avg_mass: float | None
    psi_mod: str | None
    unimod: int | None


def parse_ptmlist(text: str) -> dict[str, UniProtPtm]:
    """Parse UniProt ptmlist.txt and return ``{ptm_name: UniProtPtm}``.

    Records are separated by ``//`` lines. Each record has an ``ID`` line
    and optional cross-reference lines for PSI-MOD and UniMod.
    """
    result: dict[str, UniProtPtm] = {}

    current_id: str | None = None
    current_ac: str = ""
    current_ft: str = ""
    current_tg: str = ""
    current_mm: float | None = None
    current_ma: float | None = None
    current_psimod: str | None = None
    current_unimod: int | None = None

    for line in text.splitlines():
        if line.startswith("ID   "):
            current_id = line[5:].strip()
        elif line.startswith("AC   "):
            current_ac = line[5:].strip()
        elif line.startswith("FT   "):
            current_ft = line[5:].strip()
        elif line.startswith("TG   "):
            current_tg = line[5:].strip().rstrip(".")
        elif line.startswith("MM   "):
            try:
                current_mm = float(line[5:].strip())
            except ValueError:
                pass
        elif line.startswith("MA   "):
            try:
                current_ma = float(line[5:].strip())
            except ValueError:
                pass
        elif line.startswith("DR   PSI-MOD;"):
            # e.g. "DR   PSI-MOD; MOD:00046."
            current_psimod = line.split(";")[1].strip().rstrip(".")
        elif line.startswith("DR   Unimod;"):
            # e.g. "DR   Unimod; 35."
            try:
                current_unimod = int(line.split(";")[1].strip().rstrip("."))
            except ValueError:
                pass
        elif line.strip() == "//":
            if current_id:
                result[current_id] = UniProtPtm(
                    id=current_id,
                    ac=current_ac,
                    feature_key=current_ft,
                    target=current_tg,
                    mono_mass=current_mm,
                    avg_mass=current_ma,
                    psi_mod=current_psimod,
                    unimod=current_unimod,
                )
            current_id = None
            current_ac = ""
            current_ft = ""
            current_tg = ""
            current_mm = None
            current_ma = None
            current_psimod = None
            current_unimod = None

    return result


def get_ptm_map() -> dict[str, UniProtPtm]:
    """Fetch and return the PTM map, caching after the first call."""
    global _PTM_CACHE  # noqa: PLW0603
    if _PTM_CACHE is None:
        _PTM_CACHE = parse_ptmlist(fetch_ptmlist())
    return _PTM_CACHE
