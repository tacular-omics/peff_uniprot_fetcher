"""Convert GFF feature dicts to pefftacular annotation objects."""

from __future__ import annotations

import re

from pefftacular import ModRes, ModResPsi, ModResUnimod, Processed, VariantComplex, VariantSimple
from uniprotptmpy import PtmEntry

from peff_uniprot_fetcher._ptm import get_psimod, get_unimod, psi_mod_accession, unimod_accession

VARIANT_PATTERN = re.compile(r"([A-Z]+)\s*->\s*([A-Z]+)")
DBSNP_PATTERN = re.compile(r"dbSNP:(rs\d+)")

PROCESSED_ACCESSIONS: dict[str, tuple[str, str]] = {
    "Signal peptide": ("PEFF:0001001", "signal peptide"),
    "Transit peptide": ("PEFF:0001002", "transit peptide"),
    "Propeptide": ("PEFF:0001003", "propeptide"),
    "Chain": ("PEFF:0001004", "mature protein"),
    "Peptide": ("PEFF:0001005", "peptide"),
}

_VARIANT_FEATURES = frozenset(
    {
        "Natural variant",
        "Mutagenesis",
        "Alternative sequence",
        "Sequence conflict",
    }
)

_MODIFICATION_FEATURES = frozenset(
    {
        "Modified residue",
        "Glycosylation",
        "Lipidation",
        "Cross-link",
    }
)

_PROCESSED_FEATURES = frozenset(PROCESSED_ACCESSIONS)


def _clean_mod_name(note: str) -> str:
    """Extract a clean modification name from a GFF Note value.

    Takes the first semicolon-delimited part and strips qualifier
    suffixes like ``(by ...)``, ``(Microbial infection)``, and
    ``alternate``.  Also ensures parentheses are balanced, since PEFF
    uses ``(...)`` as item delimiters and an unclosed ``(`` in a name
    would produce unreadable output.
    """
    name = note.split(";")[0].strip()
    # Remove parenthesized qualifiers at the end.
    name = re.sub(r"\s*\(.*?\)\s*$", "", name)
    name = name.strip()
    # Safety net: strip all parens if they're unbalanced.
    if name.count("(") != name.count(")"):
        name = name.replace("(", "").replace(")", "")
    return name


def _psi_name(accession: str) -> str:
    """Return the PSI-MOD name for *accession* prefixed with ``M:``, or ``M:{accession}`` if not found."""
    info = get_psimod(accession)
    return f"M:{info.name}" if info is not None else f"M:{accession}"


def _unimod_name(accession: int) -> str:
    """Return the UniMod name for *accession* prefixed with ``U:``, or ``U:{accession}`` if not found."""
    info = get_unimod(accession)
    return f"U:{info.name}" if info is not None else f"U:{accession}"


def _resolve_modification(
    mod_name: str,
    positions: tuple[int, ...],
    ptm_map: dict[str, PtmEntry],
    only_known_mass: bool,
    mod_res_psi: list[ModResPsi],
    mod_res_unimod: list[ModResUnimod],
    mod_res: list[ModRes],
) -> None:
    """Look up *mod_name* in *ptm_map* and append resolved annotations.

    Each annotation type (PSI-MOD, UniMod, generic ModRes) is resolved
    independently so a failure in one does not block the others.  When
    no PTM match is found a generic ``ModRes`` with an empty accession
    is emitted as a fallback.
    """
    ptm = ptm_map.get(mod_name)
    psi_mod_id = psi_mod_accession(ptm) if ptm else None
    unimod_id = unimod_accession(ptm) if ptm else None

    if psi_mod_id:
        psi_mod = get_psimod(psi_mod_id)
        if psi_mod is not None:
            if not only_known_mass or psi_mod.mass_mono is not None:
                mod_res_psi.append(ModResPsi(positions=positions, accession=psi_mod_id, name=_psi_name(psi_mod_id)))

    if unimod_id:
        unimod = get_unimod(unimod_id)
        if unimod is not None:
            if not only_known_mass or unimod.delta_mono_mass is not None:
                mod_res_unimod.append(
                    ModResUnimod(positions=positions, accession=f"UNIMOD:{unimod_id}", name=_unimod_name(unimod_id))
                )

    if ptm:
        has_mass = ptm.monoisotopic_mass is not None
        if not has_mass and only_known_mass:
            return
        mod_res.append(ModRes(positions=positions, accession=ptm.id, name=ptm.name))


def features_to_annotations(
    features: list[dict],
    ptm_map: dict[str, PtmEntry],
    *,
    only_known_mass: bool = False,
) -> dict:
    """Convert GFF feature dicts to PEFF annotation tuples.

    Parameters
    ----------
    features:
        List of feature dicts as returned by :func:`_gff.parse_gff`.
    ptm_map:
        Mapping of PTM name to :class:`~uniprotptmpy.PtmEntry`.

    Returns
    -------
    dict
        Keys: ``variant_simple``, ``variant_complex``, ``mod_res_psi``,
        ``mod_res``, ``processed``.  Values are sorted tuples of the
        respective annotation objects.
    """
    variant_simple: list[VariantSimple] = []
    variant_complex: list[VariantComplex] = []
    mod_res_unimod: list[ModResUnimod] = []
    mod_res_psi: list[ModResPsi] = []
    mod_res: list[ModRes] = []
    processed: list[Processed] = []

    for feat in features:
        ftype = feat["feature"]
        start: int = feat["start"]
        end: int = feat["end"]
        attrs: dict[str, str] = feat["attributes"]
        note = attrs.get("Note", "")

        # -- Variants ---------------------------------------------------
        if ftype in _VARIANT_FEATURES:
            match = VARIANT_PATTERN.search(note)
            if match:
                original, mutant = match.group(1), match.group(2)
                # Detect tag (dbSNP ID).
                dbsnp_match = DBSNP_PATTERN.search(note)
                tag = dbsnp_match.group(1) if dbsnp_match else None

                if start == end and len(original) == 1 and len(mutant) == 1:
                    variant_simple.append(VariantSimple(position=start, new_amino_acid=mutant, tag=tag))
                else:
                    variant_complex.append(VariantComplex(start_pos=start, end_pos=end, new_sequence=mutant, tag=tag))
            elif "Missing" in note:
                variant_complex.append(VariantComplex(start_pos=start, end_pos=end, new_sequence=""))

        # -- Modified residue -------------------------------------------
        elif ftype == "Modified residue":
            mod_name = _clean_mod_name(note)
            if not mod_name:
                continue
            _resolve_modification(mod_name, (start,), ptm_map, only_known_mass, mod_res_psi, mod_res_unimod, mod_res)

        # -- Glycosylation / Lipidation ---------------------------------
        elif ftype in ("Glycosylation", "Lipidation"):
            mod_name = _clean_mod_name(note) if note else ftype
            # Try raw note first (may match PTM names that include parenthesized
            # detail, e.g. "N-linked (GlcNAc...)"), then fall back to cleaned name.
            raw_name = note.split(";")[0].strip() if note else ftype
            ptm = ptm_map.get(raw_name) or ptm_map.get(mod_name)
            if ptm:
                _resolve_modification(
                    ptm.name, (start,), ptm_map, only_known_mass, mod_res_psi, mod_res_unimod, mod_res
                )

        # -- Cross-link -------------------------------------------------
        elif ftype == "Cross-link":
            mod_name = _clean_mod_name(note) if note else "Cross-link"
            positions = (start, end) if start != end else (start,)
            mod_res.append(ModRes(positions=positions, accession="", name=mod_name))

        # -- Processed features -----------------------------------------
        elif ftype in _PROCESSED_FEATURES:
            acc, name = PROCESSED_ACCESSIONS[ftype]
            processed.append(Processed(start_pos=start, end_pos=end, accession=acc, name=name))

    return {
        "variant_simple": tuple(sorted(variant_simple, key=lambda v: v.position)),
        "variant_complex": tuple(sorted(variant_complex, key=lambda v: v.start_pos)),
        "mod_res_unimod": tuple(sorted(mod_res_unimod, key=lambda m: m.positions[0])),
        "mod_res_psi": tuple(sorted(mod_res_psi, key=lambda m: m.positions[0])),
        "mod_res": tuple(sorted(mod_res, key=lambda m: m.positions[0])),
        "processed": tuple(sorted(processed, key=lambda p: p.start_pos)),
    }
