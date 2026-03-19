"""Convert GFF feature dicts to pefftacular annotation objects."""

from __future__ import annotations

import re

from pefftacular import ModRes, ModResPsi, Processed, VariantComplex, VariantSimple

VARIANT_PATTERN = re.compile(r"([A-Z]+)\s*->\s*([A-Z]+)")
DBSNP_PATTERN = re.compile(r"dbSNP:(rs\d+)")

PROCESSED_ACCESSIONS: dict[str, tuple[str, str]] = {
    "Signal peptide": ("PEFF:0001001", "signal peptide"),
    "Transit peptide": ("PEFF:0001002", "transit peptide"),
    "Propeptide": ("PEFF:0001003", "propeptide"),
    "Chain": ("PEFF:0001004", "mature protein"),
    "Peptide": ("PEFF:0001005", "peptide"),
}

_VARIANT_FEATURES = frozenset({
    "Natural variant",
    "Mutagenesis",
    "Alternative sequence",
    "Sequence conflict",
})

_MODIFICATION_FEATURES = frozenset({
    "Modified residue",
    "Glycosylation",
    "Lipidation",
    "Cross-link",
})

_PROCESSED_FEATURES = frozenset(PROCESSED_ACCESSIONS)


def _clean_mod_name(note: str) -> str:
    """Extract a clean modification name from a GFF Note value.

    Takes the first semicolon-delimited part and strips qualifier
    suffixes like ``(by ...)``, ``(Microbial infection)``, and
    ``alternate``.
    """
    name = note.split(";")[0].strip()
    # Remove parenthesized qualifiers at the end.
    name = re.sub(r"\s*\(.*?\)\s*$", "", name)
    return name.strip()


def features_to_annotations(
    features: list[dict],
    ptm_map: dict[str, str],
) -> dict:
    """Convert GFF feature dicts to PEFF annotation tuples.

    Parameters
    ----------
    features:
        List of feature dicts as returned by :func:`_gff.parse_gff`.
    ptm_map:
        Mapping of PTM name to PSI-MOD accession from :func:`_ptm.parse_ptmlist`.

    Returns
    -------
    dict
        Keys: ``variant_simple``, ``variant_complex``, ``mod_res_psi``,
        ``mod_res``, ``processed``.  Values are sorted tuples of the
        respective annotation objects.
    """
    variant_simple: list[VariantSimple] = []
    variant_complex: list[VariantComplex] = []
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
                    variant_simple.append(
                        VariantSimple(position=start, new_amino_acid=mutant, tag=tag)
                    )
                else:
                    variant_complex.append(
                        VariantComplex(start_pos=start, end_pos=end, new_sequence=mutant, tag=tag)
                    )
            elif "Missing" in note:
                variant_complex.append(
                    VariantComplex(start_pos=start, end_pos=end, new_sequence="")
                )

        # -- Modified residue -------------------------------------------
        elif ftype == "Modified residue":
            mod_name = _clean_mod_name(note)
            if not mod_name:
                continue
            psi_acc = ptm_map.get(mod_name)
            if psi_acc:
                mod_res_psi.append(
                    ModResPsi(positions=(start,), accession=psi_acc, name=mod_name)
                )
            else:
                mod_res.append(
                    ModRes(positions=(start,), accession="", name=mod_name)
                )

        # -- Glycosylation / Lipidation ---------------------------------
        elif ftype in ("Glycosylation", "Lipidation"):
            mod_name = _clean_mod_name(note) if note else ftype
            mod_res.append(
                ModRes(positions=(start,), accession="", name=mod_name)
            )

        # -- Cross-link -------------------------------------------------
        elif ftype == "Cross-link":
            mod_name = _clean_mod_name(note) if note else "Cross-link"
            positions = (start, end) if start != end else (start,)
            mod_res.append(
                ModRes(positions=positions, accession="", name=mod_name)
            )

        # -- Processed features -----------------------------------------
        elif ftype in _PROCESSED_FEATURES:
            acc, name = PROCESSED_ACCESSIONS[ftype]
            processed.append(
                Processed(start_pos=start, end_pos=end, accession=acc, name=name)
            )

    return {
        "variant_simple": tuple(sorted(variant_simple, key=lambda v: v.position)),
        "variant_complex": tuple(sorted(variant_complex, key=lambda v: v.start_pos)),
        "mod_res_psi": tuple(sorted(mod_res_psi, key=lambda m: m.positions[0])),
        "mod_res": tuple(sorted(mod_res, key=lambda m: m.positions[0])),
        "processed": tuple(sorted(processed, key=lambda p: p.start_pos)),
    }
