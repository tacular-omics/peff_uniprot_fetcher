# History

## [Unreleased]

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
