"""PTM name to UniProt PTM entry mapping."""

from __future__ import annotations

from dataclasses import dataclass, replace

from tacular.psimod import PSIMOD_LOOKUP
from tacular.unimod import UNIMOD_LOOKUP

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
    formula:
        Elemental correction formula as listed in ptmlist.txt (e.g. ``"H1 O3 P1"``),
        or ``None`` if not listed.
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
    formula: str | None
    mono_mass: float | None
    avg_mass: float | None
    psi_mod: str | None
    unimod: int | None

    @property
    def proforma_formula(self) -> str | None:
        """ProForma 2.0 formula tag, e.g. ``[Formula:H1O3P1]``.

        Converts the space-separated UniProt CF format (``"H1 O3 P1"``) to the
        compact ProForma notation by removing spaces between element tokens.
        Returns ``None`` if no formula is available.
        """
        if self.formula is None:
            return None
        return f"Formula:{''.join(self.formula.split())}"


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
    current_cf: str | None = None
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
        elif line.startswith("CF   "):
            current_cf = line[5:].strip()
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
                    formula=current_cf,
                    mono_mass=current_mm,
                    avg_mass=current_ma,
                    psi_mod=current_psimod,
                    unimod=current_unimod,
                )
            current_id = None
            current_ac = ""
            current_ft = ""
            current_tg = ""
            current_cf = None
            current_mm = None
            current_ma = None
            current_psimod = None
            current_unimod = None

    return result


def _enrich_from_tacular(
    ptm_map: dict[str, UniProtPtm],
) -> dict[str, UniProtPtm]:
    """Fill formula / mono_mass / avg_mass from tacular for entries that
    have a PSI-MOD or UniMod cross-reference but no mass data from UniProt."""
    enriched: dict[str, UniProtPtm] = {}
    for name, ptm in ptm_map.items():
        if ptm.formula is not None or (ptm.unimod is None and ptm.psi_mod is None):
            enriched[name] = ptm
            continue

        # Prefer UniMod (more standardised masses), fall back to PSI-MOD
        info = None
        if ptm.unimod is not None:
            info = UNIMOD_LOOKUP.query_id(ptm.unimod)
        if info is None and ptm.psi_mod is not None:
            # psi_mod is stored as "MOD:00046"; the lookup wants the numeric part
            try:
                psi_mod_num = int(ptm.psi_mod.split(":")[1])
            except (IndexError, ValueError):
                psi_mod_num = None
            if psi_mod_num is not None:
                info = PSIMOD_LOOKUP.query_id(psi_mod_num)

        if info is None:
            enriched[name] = ptm
            continue

        # Convert dict_composition → UniProt "El1 El2 ..." format
        formula: str | None = None
        if info.dict_composition:
            formula = " ".join(
                f"{el}{cnt}"
                for el, cnt in sorted(info.dict_composition.items())
            )

        enriched[name] = replace(
            ptm,
            formula=formula if formula is not None else ptm.formula,
            mono_mass=info.monoisotopic_mass if ptm.mono_mass is None else ptm.mono_mass,
            avg_mass=info.average_mass if ptm.avg_mass is None else ptm.avg_mass,
        )
    return enriched


def get_ptm_map() -> dict[str, UniProtPtm]:
    """Fetch and return the PTM map, caching after the first call."""
    global _PTM_CACHE  # noqa: PLW0603
    if _PTM_CACHE is None:
        _PTM_CACHE = _enrich_from_tacular(parse_ptmlist(fetch_ptmlist()))
    return _PTM_CACHE
