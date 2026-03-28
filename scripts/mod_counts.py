"""Count modification occurrences across all entries in a PEFF file."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from pefftacular import PeffReader

from peff_uniprot_fetcher._ptm import get_ptm_map, psi_mod_accession, unimod_accession


def _build_mass_lookups(ptm_map: dict) -> tuple[dict[str, bool], dict[str, bool]]:
    """Return (psi_accession -> has_mass, unimod_accession -> has_mass) dicts."""
    psi = {}
    unimod = {}
    for v in ptm_map.values():
        psi_acc = psi_mod_accession(v)
        uni_acc = unimod_accession(v)
        if psi_acc:
            psi[psi_acc] = v.monoisotopic_mass is not None
        if uni_acc is not None:
            unimod[f"UNIMOD:{uni_acc}"] = v.monoisotopic_mass is not None
    return psi, unimod


def count_mods(peff_path: Path) -> dict[str, dict]:
    """Return a dict mapping mod name -> {total, mod_res_psi, mod_res_unimod, mod_res, accessions}."""
    counts: dict[str, Counter] = defaultdict(Counter)
    accessions: dict[str, set] = defaultdict(set)

    with PeffReader(peff_path) as reader:
        for entry in reader:
            for m in entry.mod_res_psi:
                counts[m.name]["mod_res_psi"] += 1
                accessions[m.name].add(("psi", m.accession))
            for m in entry.mod_res_unimod:
                counts[m.name]["mod_res_unimod"] += 1
                accessions[m.name].add(("unimod", m.accession))
            for m in entry.mod_res:
                counts[m.name]["mod_res"] += 1
                accessions[m.name].add(("mod", m.accession))

    return {
        name: {
            "mod_res_psi": c["mod_res_psi"],
            "mod_res_unimod": c["mod_res_unimod"],
            "mod_res": c["mod_res"],
            "total": c["mod_res_psi"] + c["mod_res_unimod"] + c["mod_res"],
            "accessions": accessions[name],
        }
        for name, c in counts.items()
    }


def write_tsv(path: Path, mods: dict[str, dict]) -> None:
    ptm_map = get_ptm_map()
    psi_mass, unimod_mass = _build_mass_lookups(ptm_map)

    rows = []
    for name, data in mods.items():
        has_mass = False
        for kind, acc in data["accessions"]:
            if kind == "psi" and psi_mass.get(acc):
                has_mass = True
                break
            if kind == "unimod" and unimod_mass.get(acc):
                has_mass = True
                break
            if kind == "mod" and acc.startswith("Formula:"):
                has_mass = True
                break

        rows.append({
            "name": name,
            "total": data["total"],
            "mod_res_psi": data["mod_res_psi"],
            "mod_res_unimod": data["mod_res_unimod"],
            "mod_res": data["mod_res"],
            "has_mass": "true" if has_mass else "false",
        })

    rows.sort(key=lambda r: r["total"], reverse=True)

    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["name", "total", "mod_res_psi", "mod_res_unimod", "mod_res", "has_mass"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)

    print(f"Written {len(rows)} modifications to {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count modification occurrences in a PEFF file.")
    parser.add_argument("input", metavar="INPUT", help="Path to input PEFF file.")
    parser.add_argument("output", metavar="OUTPUT", help="Path for output TSV file.")
    args = parser.parse_args()

    mods = count_mods(Path(args.input))

    total = sum(d["total"] for d in mods.values())
    psi_total = sum(d["mod_res_psi"] for d in mods.values())
    unimod_total = sum(d["mod_res_unimod"] for d in mods.values())
    mod_total = sum(d["mod_res"] for d in mods.values())
    print(f"Total modification annotations: {total}")
    print(f"  ModResPsi:    {psi_total}  ({sum(1 for d in mods.values() if d['mod_res_psi'])} distinct)")
    print(f"  ModResUnimod: {unimod_total}  ({sum(1 for d in mods.values() if d['mod_res_unimod'])} distinct)")
    print(f"  ModRes:       {mod_total}  ({sum(1 for d in mods.values() if d['mod_res'])} distinct)")
    print()

    write_tsv(Path(args.output), mods)
