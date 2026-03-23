"""Generate a PEFF file for the reviewed human proteome and report modification statistics.

Statistics are collected on actual annotation instances (one per annotated site, not per
unique PTM name), broken down by feature type (MOD_RES / CROSSLNK / LIPID / CARBOHYD)
and ontology mapping (PSI-MOD / UniMod / both / formula-only / none).
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

from pefftacular import write_peff

from peff_uniprot_fetcher._annotations import _clean_mod_name, features_to_annotations
from peff_uniprot_fetcher._builder import build_entry, build_header
from peff_uniprot_fetcher._client import stream_search
from peff_uniprot_fetcher._fasta import parse_fasta
from peff_uniprot_fetcher._gff import parse_gff
from peff_uniprot_fetcher._ptm import UniProtPtm, get_ptm_map

log = logging.getLogger(__name__)

HUMAN_QUERY = "reviewed:true AND organism_id:9606"

_MOD_FEATURE_KEYS: dict[str, str] = {
    "Modified residue": "MOD_RES",
    "Glycosylation": "CARBOHYD",
    "Lipidation": "LIPID",
    "Cross-link": "CROSSLNK",
}

_VARIANT_TYPES = {"Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict"}
_PROCESSED_TYPES = {"Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"}

_MAPPING_TYPES = ("psi_only", "unimod_only", "both", "custom", "none")


def _mapping_type(ptm: UniProtPtm | None) -> str:
    if ptm is None:
        return "none"
    has_psi = ptm.psi_mod is not None
    has_uni = ptm.unimod is not None
    if has_psi and has_uni:
        return "both"
    if has_psi:
        return "psi_only"
    if has_uni:
        return "unimod_only"
    if ptm.proforma_formula is not None:
        return "custom"
    return "none"


def collect_mod_stats(
    all_features: dict[str, list[dict]],
    ptm_map: dict[str, UniProtPtm],
) -> tuple[dict[str, Counter], Counter, Counter]:
    """Walk all protein features and count modification instances.

    Returns
    -------
    per_fk : dict[feature_key, Counter]
        Counts keyed by mapping type + ``"has_mass"``.
    mod_name_counts : Counter
        How many times each PTM name was observed across all proteins.
    unmapped_counts : Counter
        PTM names that had no entry in *ptm_map*, with occurrence counts.
    """
    per_fk: dict[str, Counter] = defaultdict(Counter)
    mod_name_counts: Counter = Counter()
    unmapped_counts: Counter = Counter()

    for features in all_features.values():
        for feat in features:
            ftype = feat["feature"]
            if ftype not in _MOD_FEATURE_KEYS:
                continue
            fk = _MOD_FEATURE_KEYS[ftype]
            note = feat["attributes"].get("Note", "")
            mod_name = _clean_mod_name(note) if note else ftype

            ptm = ptm_map.get(mod_name)
            mtype = _mapping_type(ptm)
            has_mass = ptm is not None and ptm.mono_mass is not None

            per_fk[fk][mtype] += 1
            if has_mass:
                per_fk[fk]["has_mass"] += 1

            mod_name_counts[mod_name] += 1
            if ptm is None:
                unmapped_counts[mod_name] += 1

    return per_fk, mod_name_counts, unmapped_counts


def print_stats(
    per_fk: dict[str, Counter],
    mod_name_counts: Counter,
    unmapped_counts: Counter,
    top_n: int = 20,
) -> None:
    fk_order = ["MOD_RES", "CARBOHYD", "CROSSLNK", "LIPID"]
    all_fks = fk_order + [fk for fk in sorted(per_fk) if fk not in fk_order]

    print("=" * 60)
    print("Modification instance statistics (annotated proteome)")
    print("=" * 60)

    grand_total = sum(sum(c.values()) - c["has_mass"] for c in per_fk.values())
    grand_mass  = sum(c["has_mass"] for c in per_fk.values())

    print(f"\nTotal modification instances: {grand_total + grand_mass:,}")
    print(f"  With monoisotopic mass:      {grand_mass:,}  ({grand_mass / max(grand_total + grand_mass, 1):.1%})")
    print()

    header = f"{'Feature':10}  {'Total':>8}  {'has_mass':>9}  {'PSI-only':>9}  {'Uni-only':>9}  {'Both':>6}  {'Custom':>7}  {'None':>6}"
    print(header)
    print("-" * len(header))

    for fk in all_fks:
        if fk not in per_fk:
            continue
        c = per_fk[fk]
        total = sum(v for k, v in c.items() if k != "has_mass")
        print(
            f"{fk:10}  {total:>8,}  {c['has_mass']:>9,}  "
            f"{c['psi_only']:>9,}  {c['unimod_only']:>9,}  "
            f"{c['both']:>6,}  {c['custom']:>7,}  {c['none']:>6,}"
        )

    print()
    print(f"Top {top_n} most common modification names (across all proteins):")
    for name, count in mod_name_counts.most_common(top_n):
        tag = " [unmapped]" if name in unmapped_counts else ""
        print(f"  {count:>8,}  {name}{tag}")

    if unmapped_counts:
        print()
        print(f"Top {top_n} unmapped modification names (no entry in ptm_map):")
        for name, count in unmapped_counts.most_common(top_n):
            print(f"  {count:>8,}  {name}")


def run(
    output: Path,
    query: str,
    include_variants: bool,
    include_modifications: bool,
    include_glycosylation: bool,
    include_lipidation: bool,
    include_crosslinks: bool,
    include_processed: bool,
    exclusive_mod_lists: bool,
    top_n: int,
) -> None:
    log.info("Fetching FASTA for query: %s", query)
    fasta_entries = parse_fasta(stream_search(query, fmt="fasta"))  # type: ignore[arg-type]
    log.info("Parsed %d sequences", len(fasta_entries))

    log.info("Fetching GFF annotations (this may take a while)...")
    all_features = parse_gff(stream_search(query, fmt="gff"))  # type: ignore[arg-type]
    log.info("Parsed GFF for %d accessions", len(all_features))

    ptm_map = get_ptm_map()
    log.info("PTM map loaded (%d entries)", len(ptm_map))

    # Collect stats from raw features before building PEFF entries.
    log.info("Collecting modification statistics...")
    per_fk, mod_name_counts, unmapped_counts = collect_mod_stats(all_features, ptm_map)

    # Build PEFF entries.
    log.info("Building PEFF entries...")
    _variant_types  = {"Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict"}
    _processed_types = {"Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"}

    entries = []
    for fasta_entry in fasta_entries:
        raw = all_features.get(fasta_entry.accession, [])
        filtered: list[dict] = []
        for feat in raw:
            ft = feat["feature"]
            if ft in _variant_types and include_variants:
                filtered.append(feat)
            elif ft == "Modified residue" and include_modifications:
                filtered.append(feat)
            elif ft == "Glycosylation" and include_glycosylation:
                filtered.append(feat)
            elif ft == "Lipidation" and include_lipidation:
                filtered.append(feat)
            elif ft == "Cross-link" and include_crosslinks:
                filtered.append(feat)
            elif ft in _processed_types and include_processed:
                filtered.append(feat)

        annotations = features_to_annotations(filtered, ptm_map, exclusive_mod_lists=exclusive_mod_lists)
        entries.append(build_entry(fasta_entry, annotations))

    header = build_header(entries)
    log.info("Writing PEFF to %s...", output)
    write_peff(header, entries, output)
    log.info("Done. Wrote %d entries.", len(entries))

    print_stats(per_fk, mod_name_counts, unmapped_counts, top_n=top_n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a PEFF file for the human proteome and report modification stats."
    )
    parser.add_argument(
        "output", nargs="?", default="human_proteome.peff",
        metavar="OUTPUT",
        help="Path to write the PEFF file (default: human_proteome.peff).",
    )
    parser.add_argument(
        "--query", default=HUMAN_QUERY, metavar="QUERY",
        help=f"UniProt search query (default: {HUMAN_QUERY!r}).",
    )
    parser.add_argument("--no-variants",           dest="variants",      action="store_false", default=True)
    parser.add_argument("--no-modifications",      dest="modifications", action="store_false", default=True)
    parser.add_argument("--include-glycosylation", dest="glycosylation", action="store_true",  default=False)
    parser.add_argument("--include-lipidation",    dest="lipidation",    action="store_true",  default=False)
    parser.add_argument("--include-crosslinks",    dest="crosslinks",    action="store_true",  default=False)
    parser.add_argument("--no-processed",          dest="processed",     action="store_false", default=True)
    parser.add_argument(
        "--exclusive-mod-lists", action="store_true", default=False,
        help="Assign each modification to exactly one list (PSI > UniMod > Custom).",
    )
    parser.add_argument(
        "--top-n", type=int, default=20, metavar="N",
        help="Number of top modification names to display (default: 20).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s", stream=sys.stderr)

    run(
        output=Path(args.output),
        query=args.query,
        include_variants=args.variants,
        include_modifications=args.modifications,
        include_glycosylation=args.glycosylation,
        include_lipidation=args.lipidation,
        include_crosslinks=args.crosslinks,
        include_processed=args.processed,
        exclusive_mod_lists=args.exclusive_mod_lists,
        top_n=args.top_n,
    )
