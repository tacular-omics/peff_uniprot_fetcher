# peff_uniprot_fetcher — Claude Code Guide

## Project overview

peff_uniprot_fetcher turns UniProt accessions, taxonomy IDs, search queries, or a local
UniProt FASTA into annotated PEFF (PSI Extended FASTA Format) files for proteomics search
engines. It downloads FASTA + GFF3 from the UniProt REST API
(`https://rest.uniprot.org/uniprotkb`), maps GFF features (variants, modified residues,
glycosylation, lipidation, cross-links, processed forms) to PEFF annotations, and writes
them with `pefftacular`. It ships three CLIs and a Pyodide web app on GitHub Pages
(<https://tacular-omics.github.io/peff_uniprot_fetcher/>, source in `docs/`).

Place in the tacular-omics graph: **downstream only**. It depends on PyPI releases of
`pefftacular` (PEFF models + writer), `uniprotptmpy` (UniProt ptmlist), `psimodpy`
(PSI-MOD) and `unimodpy` (UNIMOD), plus `httpx`. It does **not** depend on `peptacular`
or `tacular`. Nothing in the org depends on it.

## Commands

All verified with `just` / `uv run` in a fresh clone.

```bash
just install          # uv sync (creates .venv + uv.lock; uv.lock is gitignored)
just lint             # ruff check src
just format           # ruff isort fix + ruff format src (modifies files)
just ty               # ty check src
just test             # pytest tests (52 tests, offline, <1 s)
just check            # lint + ty + test
just build            # uv build
just clean            # rm -rf dist
just download-ecoli   # NETWORK: raw E. coli K-12 FASTA+GFF -> data/ecoli/
just fetch-ecoli      # NETWORK: PEFF for E. coli K-12 -> data/ecoli/ecoli.peff
just fasta-to-peff-ecoli [fasta=...]   # NETWORK: annotate the downloaded FASTA
```

CLIs (console scripts from `[project.scripts]`), all hit UniProt:

```bash
uv run fetch-peff OUT.peff (--organism-id TAXID | --query Q | --accessions ACC...) [--unreviewed] [annotation flags]
uv run fasta-to-peff IN.fasta OUT.peff [annotation flags]
uv run download-uniprot (--organism-id TAXID | --query Q | --accession ACC) [-o DIR] [--unreviewed] [--formats fasta gff]
```

Annotation flags (shared): `--no-variants --no-modifications --no-processed
--glycosylation --lipidation --crosslinks --only-known-mass`.

Scripts (not packaged): `scripts/human_proteome_peff.py` (human PEFF + mod stats, network),
`scripts/mod_counts.py` (count mods in an existing PEFF), `scripts/ptm_stats.py`
(ptmlist statistics).

## Architecture

```
src/peff_uniprot_fetcher/
  __init__.py       public API: fetch_peff / fasta_to_peff (+ *_to_file / *_file),
                    GFF batching by accession, feature filtering by AnnotationConfig
  _config.py        AnnotationConfig frozen dataclass (which feature types to emit)
  _client.py        UniProt REST client: fetch_entry, fetch_entries (/search), stream_search (/stream);
                    httpx natively, synchronous XHR under Pyodide (sys.platform == "emscripten")
  _fasta.py         UniProtFastaEntry dataclass + parse_fasta (sp|ACC|NAME ... OS= OX= GN= PE= SV=)
  _gff.py           parse_gff -> {accession: [{"feature","start","end","attributes"}]};
                    keeps only variant / modification / processed feature types
  _ptm.py           get_ptm_map(): uniprotptmpy entries by name, masses/formulas filled
                    from unimodpy (preferred) then psimodpy; lazily cached module globals
  _annotations.py   features_to_annotations(): GFF dicts -> pefftacular VariantSimple,
                    VariantComplex, ModResPsi, ModResUnimod, ModRes, Processed
  _builder.py       build_entry (SequenceEntry from FASTA + annotations), build_header
                    (one DatabaseHeader per sp/tr prefix, PEFF 1.0)
  _web.py           generate_peff_string(organism_id, reviewed=True) for the Pyodide worker
  _cli.py           argparse entry points fetch_peff_cli / fasta_to_peff_cli / download_uniprot_cli
docs/               static GitHub Pages app (index.html, app.js, worker.js) + a committed wheel
                    that worker.js installs via micropip
```

Data flow: FASTA text -> `parse_fasta`; GFF text -> `parse_gff`; `_build_entries` filters
features per `AnnotationConfig`, calls `features_to_annotations` with `get_ptm_map()`,
then `build_entry` per protein and `build_header`; `write_peff` (pefftacular) serialises.

- `fetch_peff(query=...)`: FASTA and GFF each via one `/stream` call.
- `fetch_peff(accessions=...)`: FASTA via `/search` (`fetch_entries`), GFF via
  `_fetch_gff_per_accession`, which drops non-UniProt accessions (regex, warns with
  examples). Both batch `accession:X OR ...` queries with `_client.accession_batches`
  (under 1800 characters and at most 500 accessions per request).
- `fasta_to_peff(path)`: sequences from the file, GFF batched per accession as above.

## Public API

From `peff_uniprot_fetcher.__all__`:

- `fetch_peff(accessions=None, query=None, cfg=None, **kwargs) -> (FileHeader, list[SequenceEntry])`
  exactly one of `accessions` / `query`, else `ValueError`.
- `fetch_peff_to_file(output, accessions=None, query=None, cfg=None, **kwargs) -> None`
- `fasta_to_peff(fasta, cfg=None, **kwargs) -> (FileHeader, list[SequenceEntry])`
- `fasta_to_peff_file(fasta, output, cfg=None, **kwargs) -> None`
- `AnnotationConfig(include_variants=True, include_modifications=True,
  include_glycosylation=False, include_lipidation=False, include_crosslinks=False,
  include_processed=True, only_known_mass=False)`
- Re-exported from pefftacular: `FileHeader`, `SequenceEntry`, `write_peff(header, entries, dest)`.

`**kwargs` are `AnnotationConfig` fields and are used only when `cfg is None`.

## Conventions

- Python >= 3.12, `from __future__ import annotations`, `X | None` not `Optional`,
  frozen/slotted dataclasses for config. ruff line length 120, rules E,W,F,I,B,UP.
- Docstrings: short; the modules use numpydoc-style `Parameters` / `Returns` sections.
  Match the surrounding module.
- Logging via `logging.getLogger(__name__)`; the CLIs call `logging.basicConfig(INFO)`
  to stderr. No `print` in `src/`.
- Errors: `ValueError` for bad arguments; HTTP errors propagate from
  `httpx.Response.raise_for_status()` (or `RuntimeError` under Pyodide).
- Tests in `tests/`, one file per module (`test_fasta.py`, `test_gff.py`,
  `test_annotations.py`, `test_ptm.py`, `test_cli.py`, `test_basic.py`). They are offline:
  inline FASTA/GFF strings and hand-built `PtmEntry` maps (`test_annotations.py`);
  `test_ptm.py` also loads the real `get_ptm_map()` from the installed ontology packages
  (psimodpy / unimodpy lookups for names and masses always use the real databases). Never add a test that calls UniProt.
- Lint and ty run on `src` only (CI: `just lint`, `uv run ty check src`, `pytest tests`
  on Python 3.12 and 3.13).

## Gotchas

- **No committed lock.** `uv.lock` is gitignored and CI runs `uv sync` unlocked, so every
  CI run and install picks the newest `pefftacular` / `psimodpy` / `unimodpy` /
  `uniprotptmpy`. Sibling pins are floors only (`pefftacular` has none). A breaking
  release in any of them breaks this package with no code change here.
- **`M:` / `U:` name prefixes.** `ModResPsi` names are written `M:<PSI-MOD name>` and
  `ModResUnimod` names `U:<UNIMOD name>` (e.g. `(2|MOD:00046|M:O-phospho-L-serine)`).
  Tests assert this. The README example shows the name without the prefix.
- A modified residue whose cleaned name is not in the ptmlist emits **nothing**
  (no fallback `ModRes`), despite what the `_resolve_modification` docstring says.
  Cross-links always emit a `ModRes` with an empty accession.
- Glycosylation / lipidation match the raw Note first (`N-linked (GlcNAc...) asparagine`)
  and then the qualifier-stripped name. `_clean_mod_name` strips a trailing `(...)` and
  removes unbalanced parens because PEFF uses `(...)` as delimiters.
- `_UNIPROT_ACCESSION_RE` accepts isoform accessions (`P04637-2`) and they are queried,
  but UniProt returns no GFF features for non-canonical isoforms, so they get no
  annotations (a warning names them). `parse_fasta` raises
  `ValueError` on any header without two `|`.
- `parse_gff` splits attributes on `;` only at parenthesis depth 0, after URL-decoding.
- `get_ptm_map()` and the psimod/unimod databases are cached in module globals; the first
  call loads three ontologies (about 0.3 s).
- `_client.py` must stay importable under Pyodide: `httpx` is imported only when
  `sys.platform != "emscripten"`. Keep the whole call path synchronous; the web worker
  relies on synchronous XHR.
- `docs/worker.js` installs the committed wheel `docs/peff_uniprot_fetcher-0.1.0-py3-none-any.whl`
  by filename. The web app does not pick up source changes until that wheel is rebuilt
  and the filename in `worker.js` updated.
- Large proteomes (human with GFF) are hundreds of MB and several minutes; the
  `/stream` call is one request with a 60 s timeout.
- `just argc FASTA [OUT]` runs `fasta-to-peff --only-known-mass` on a local FASTA (network).

## Releasing

Only the tacular-omics overseer bumps versions or publishes. See `just --list`.
Version source: `version` in `pyproject.toml` (static, not hatch-dynamic). Changelog:
`HISTORY.md`. A GitHub release (`release: published`) triggers
`.github/workflows/python-publish.yml`, which builds and uploads to PyPI
(`peff-uniprot-fetcher`).

## Workspace note

This repo is **not** part of the tacular-omics uv workspace
(`~/Repos/tacular-omics/packages/`). It resolves its tacular-omics dependencies from
PyPI, so run `uv run` / `just` inside this repo's own checkout. To test against an
unreleased sibling, add it explicitly (`uv run --with <path-to-sibling-checkout> pytest`).
