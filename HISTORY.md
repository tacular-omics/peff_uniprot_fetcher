# History

## [Unreleased]

### Changed

* Requires the 1.x releases of the PEFF and modification libraries:
  `pefftacular>=1.0,<2`, `psimodpy>=1.0,<2`, `unimodpy>=1.0,<2`, `uniprotptmpy>=1.0,<2`
  (was the 0.x series). No output change.
* PSI-MOD formulas are read from `dict_composition` instead of the deprecated
  `dict_diff_formula` alias (removed in psimodpy 2.0).

## 0.2.0 (2026-09-23)

### Changed

* **Output format change (PEFF 1.0 compliance).** `\ModResPsi` and `\ModResUnimod`
  names are now the bare ontology `name:`, as PEFF 1.0 sections 3.3.10/3.3.11 require.
  The `M:` / `U:` prefixes are gone: `(2|MOD:00046|M:O-phospho-L-serine)` is now
  `(2|MOD:00046|O-phospho-L-serine)` and `(2|UNIMOD:21|U:Phospho)` is now
  `(2|UNIMOD:21|Phospho)`.
* **Output format change.** A `\ModRes` entry (UniProt ptmlist id) is now written only
  when the modification maps to neither PSI-MOD nor UNIMOD (PEFF 1.0 section 3.3.12).
  Previously every ptmlist match also got a `\ModRes`, duplicating the CV entries.
  Cross-links are unchanged.

### Fixed

* Isoform accessions (e.g. `P04637-2`) are no longer rejected as non-UniProt. They are
  queried like any other accession; because UniProt publishes GFF features only for
  canonical sequences, a warning now names every isoform that got no annotations.
  Accessions that are still rejected (e.g. contaminants) are named in the warning too.
* `fetch_peff(accessions=...)` now pages the FASTA request: accessions are sent in
  batches of at most 500 and under UniProt's query-length limit. Previously lists over
  500 were silently truncated.
* `only_known_mass` now checks the PSI-MOD delta mass (`diff_mono`), consistent with the
  UNIMOD delta mass, instead of the full residue mass (`mass_mono`).
* Added `[project.urls]` (Homepage, Repository, Issues, Changelog) to `pyproject.toml`.
* The `just argc` recipe takes the FASTA path as an argument instead of a hard-coded
  personal path.

## 0.1.1 (2026-09-23)

* Archived on Zenodo (`CITATION.cff`, `.zenodo.json`).
* Capped sibling dependency pins to the org's floor-plus-major-cap policy:
  `pefftacular>=0.4,<0.5`, `psimodpy>=0.2,<0.3`, `unimodpy>=0.2,<0.3`,
  `uniprotptmpy>=0.2,<0.3` (previously `>=0.1.2`/unpinned and uncapped).
  No code changes were needed; the existing test suite passes unchanged
  against the new floors.

## 0.1.0 (2026-03-18)

* First release on PyPI.
