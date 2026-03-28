"""Print statistics on UniProt's ptmlist.txt entries and write subset TSV files."""

import csv
import sys
from collections import Counter
from pathlib import Path

from uniprotptmpy import PtmEntry

from peff_uniprot_fetcher._ptm import get_ptm_map, psi_mod_accession, unimod_accession

FIELDS = [
    "name", "id", "feature_type", "target", "correction_formula",
    "proforma_formula", "monoisotopic_mass", "average_mass", "psi_mod", "unimod",
]


def write_tsv(path: Path, entries: list[PtmEntry]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for p in sorted(entries, key=lambda x: x.name):
            w.writerow({
                "name": p.name,
                "id": p.id,
                "feature_type": p.feature_type,
                "target": p.target,
                "correction_formula": p.correction_formula or "",
                "proforma_formula": p.proforma_formula or "",
                "monoisotopic_mass": "" if p.monoisotopic_mass is None else p.monoisotopic_mass,
                "average_mass": "" if p.average_mass is None else p.average_mass,
                "psi_mod": psi_mod_accession(p) or "",
                "unimod": unimod_accession(p) or "",
            })
    print(f"  {path}  ({len(entries)} entries)")


def write_subsets(ptms: dict[str, PtmEntry], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = list(ptms.values())

    subsets = {
        "all":          entries,
        "psimod_only":  [p for p in entries if psi_mod_accession(p) is not None and unimod_accession(p) is None],
        "unimod_only":  [p for p in entries if unimod_accession(p) is not None and psi_mod_accession(p) is None],
        "both":         [p for p in entries if psi_mod_accession(p) is not None and unimod_accession(p) is not None],
        "neither":      [p for p in entries if psi_mod_accession(p) is None and unimod_accession(p) is None],
        "MOD_RES":      [p for p in entries if p.feature_type == "MOD_RES"],
        "CARBOHYD":     [p for p in entries if p.feature_type == "CARBOHYD"],
        "CROSSLNK":     [p for p in entries if p.feature_type == "CROSSLNK"],
        "LIPID":        [p for p in entries if p.feature_type == "LIPID"],
    }

    print("Writing TSV files:")
    for name, subset in subsets.items():
        write_tsv(out_dir / f"ptm_{name}.tsv", subset)


def stats(ptms: dict[str, PtmEntry]) -> None:
    total = len(ptms)
    entries = list(ptms.values())

    has_psimod  = sum(1 for p in entries if psi_mod_accession(p) is not None)
    has_unimod  = sum(1 for p in entries if unimod_accession(p) is not None)
    has_both    = sum(1 for p in entries if psi_mod_accession(p) is not None and unimod_accession(p) is not None)
    has_neither = sum(1 for p in entries if psi_mod_accession(p) is None and unimod_accession(p) is None)
    has_mass    = sum(1 for p in entries if p.monoisotopic_mass is not None)
    ft_counts   = Counter(p.feature_type for p in entries)

    print(f"Total PTM entries:          {total}")
    print()
    print("Ontology cross-references")
    print(f"  PSI-MOD only:             {has_psimod - has_both}")
    print(f"  UniMod only:              {has_unimod - has_both}")
    print(f"  Both PSI-MOD and UniMod:  {has_both}")
    print(f"  Neither:                  {has_neither}")
    print(f"  Has PSI-MOD (total):      {has_psimod}  ({has_psimod/total:.1%})")
    print(f"  Has UniMod  (total):      {has_unimod}  ({has_unimod/total:.1%})")
    print()
    print("Monoisotopic mass available")
    print(f"  With mass:                {has_mass}  ({has_mass/total:.1%})")
    print(f"  Without mass:             {total - has_mass}")
    print()
    print("Feature key breakdown")
    for ft, count in ft_counts.most_common():
        print(f"  {ft:<20} {count}")
    print()

    masses = [p.monoisotopic_mass for p in entries if p.monoisotopic_mass is not None]
    if masses:
        masses.sort()
        n = len(masses)
        print("Monoisotopic mass distribution (Da)")
        print(f"  Min:    {masses[0]:.4f}")
        print(f"  Median: {masses[n // 2]:.4f}")
        print(f"  Max:    {masses[-1]:.4f}")
        print(f"  >0 (additions):  {sum(1 for m in masses if m > 0)}")
        print(f"  <0 (removals):   {sum(1 for m in masses if m < 0)}")
        print()


if __name__ == "__main__":
    import argparse
    import logging

    parser = argparse.ArgumentParser(description="Print UniProt PTM list statistics and write subset TSV files.")
    parser.add_argument("--out-dir", "-o", default="data/ptm", metavar="DIR",
                        help="Directory to write TSV subsets into (default: data/ptm).")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)

    ptms = get_ptm_map()
    print()
    stats(ptms)
    write_subsets(ptms, Path(args.out_dir))
