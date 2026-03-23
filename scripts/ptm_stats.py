"""Print statistics on UniProt's ptmlist.txt entries."""

import sys
from collections import Counter

from peff_uniprot_fetcher._ptm import UniProtPtm, get_ptm_map


def stats(ptms: dict[str, UniProtPtm]) -> None:
    total = len(ptms)
    entries = list(ptms.values())

    has_psimod   = sum(1 for p in entries if p.psi_mod is not None)
    has_unimod   = sum(1 for p in entries if p.unimod is not None)
    has_both     = sum(1 for p in entries if p.psi_mod is not None and p.unimod is not None)
    has_neither  = sum(1 for p in entries if p.psi_mod is None and p.unimod is None)
    has_mass     = sum(1 for p in entries if p.mono_mass is not None)

    ft_counts    = Counter(p.feature_key for p in entries)

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

    # Mass distribution for entries that have masses
    masses = [p.mono_mass for p in entries if p.mono_mass is not None]
    if masses:
        masses.sort()
        n = len(masses)
        print("Monoisotopic mass distribution (Da)")
        print(f"  Min:    {masses[0]:.4f}")
        print(f"  Median: {masses[n // 2]:.4f}")
        print(f"  Max:    {masses[-1]:.4f}")
        positives = sum(1 for m in masses if m > 0)
        negatives = sum(1 for m in masses if m < 0)
        print(f"  >0 (additions):  {positives}")
        print(f"  <0 (removals):   {negatives}")


if __name__ == "__main__":
    if "--help" in sys.argv or "-h" in sys.argv:
        print("Usage: uv run python scripts/ptm_stats.py")
        print("Fetches ptmlist.txt from UniProt and prints statistics.")
        sys.exit(0)

    import logging
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)

    print("Fetching ptmlist.txt...", file=sys.stderr)
    ptms = get_ptm_map()
    print()
    stats(ptms)
