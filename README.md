# peff_uniprot_fetcher

[![Python Package](https://github.com/tacular-omics/peff_uniprot_fetcher/actions/workflows/python-package.yml/badge.svg)](https://github.com/tacular-omics/peff_uniprot_fetcher/actions/workflows/python-package.yml)
[![PyPI](https://img.shields.io/pypi/v/peff_uniprot_fetcher)](https://pypi.org/project/peff_uniprot_fetcher/)
[![License](https://img.shields.io/github/license/tacular-omics/peff_uniprot_fetcher)](https://github.com/tacular-omics/peff_uniprot_fetcher/blob/main/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/peff_uniprot_fetcher)](https://pypi.org/project/peff_uniprot_fetcher/)

Turns UniProt accessions, taxonomy IDs, or search queries into annotated [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format) files, ready for proteomics search engines that support the format. It fetches sequences and GFF feature data from the UniProt REST API and writes them out as PEFF using [pefftacular](https://github.com/tacular-omics/pefftacular), so variants, PTMs, and processed forms end up as structured annotations instead of something you have to reconstruct from raw UniProt files yourself.

> **Try it in the browser — no install required:** <https://tacular-omics.github.io/peff_uniprot_fetcher/>

## Highlights

- **No install needed to try it** — a browser-based web app (below) generates PEFF files entirely client-side.
- **Fetch by taxonomy ID, accession list, or a raw UniProt query** — whatever fits your workflow.
- **Annotations resolved to real ontology entries**, not raw UniProt text — PTMs are matched against [psimodpy](https://pypi.org/project/psimodpy/), [unimodpy](https://pypi.org/project/unimodpy/), and [uniprotptmpy](https://pypi.org/project/uniprotptmpy/) for canonical names and masses.
- **CLI, Python API, or local FASTA conversion** — fetch straight from UniProt, or annotate a FASTA file you already have.
- **Selective annotations** — turn variants, modifications, or processed forms on or off, including opt-in glycosylation, lipidation, and cross-link support.

## Web app

A static, zero-backend web app is hosted on GitHub Pages. Type an NCBI taxonomy ID (e.g. `83333` for *E. coli* K-12, `9606` for human), click **Generate PEFF**, and the browser fetches the UniProt data, builds an annotated PEFF file, and hands you a download link — all without touching a server.

It runs the same `peff_uniprot_fetcher` Python package you'd use from the CLI, compiled to WebAssembly via [Pyodide](https://pyodide.org) inside a Web Worker. UniProt is called directly from the browser, so the page is fully client-side.

- **Open it:** <https://tacular-omics.github.io/peff_uniprot_fetcher/>
- **When to use the CLI instead:** large proteomes (human, mouse, plants) pull hundreds of MB of GFF and take several minutes in-tab — use the CLI or Python API below for those. The web app is best for small/medium organisms and quick one-offs.
- **Source:** the static bundle lives in [`docs/`](https://github.com/tacular-omics/peff_uniprot_fetcher/tree/main/docs) and loads the project wheel via `micropip`; see `docs/worker.js` for the boot sequence.

## Installation

Only needed if you want the CLI or Python API — skip this if you're using the web app above.

```bash
# From PyPI (recommended)
pip install peff_uniprot_fetcher
# or, with uv
uv pip install peff_uniprot_fetcher
```

```bash
# From source
git clone https://github.com/tacular-omics/peff_uniprot_fetcher
cd peff_uniprot_fetcher
just install
```

## CLI Usage

### Fetch PEFF by organism

```bash
# Human Swiss-Prot proteome (reviewed only)
fetch-peff human.peff --organism-id 9606

# E. coli K-12
fetch-peff ecoli.peff --organism-id 83333

# Include unreviewed (TrEMBL) entries
fetch-peff human_full.peff --organism-id 9606 --unreviewed

# Custom UniProt query
fetch-peff kinases.peff --query "organism_id:9606 AND keyword:KW-0418"

# Specific accessions
fetch-peff selected.peff --accessions P12345 Q99999 O75807

# Sequences only, no annotations
fetch-peff seqs.peff --organism-id 9606 --no-variants --no-modifications --no-processed
```

### Convert a local FASTA to PEFF

Sequences come from the local file; GFF annotations are fetched from UniProt per accession.

```bash
fasta-to-peff input.fasta output.peff
```

### Download raw UniProt files

Download FASTA and/or GFF files for local inspection.

```bash
# Single accession
download-uniprot --accession P04637

# Full organism (both formats)
download-uniprot --organism-id 9606 --output-dir data/human

# GFF only
download-uniprot --organism-id 9606 --formats gff --output-dir data/human
```

### Annotation flags

All `fetch-peff` and `fasta-to-peff` commands accept:

| Flag | Default | Effect |
|---|---|---|
| `--no-variants` | variants on | Exclude sequence variants (`VariantSimple`, `VariantComplex`) |
| `--no-modifications` | modifications on | Exclude PTMs (`ModResPsi`, `ModResUnimod`, `ModRes`) |
| `--no-processed` | processed on | Exclude processed forms (`Signal peptide`, `Chain`, etc.) |
| `--glycosylation` | off | Include glycosylation annotations |
| `--lipidation` | off | Include lipidation annotations |
| `--crosslinks` | off | Include cross-link annotations |
| `--only-known-mass` | off | Only include modifications with a known monoisotopic mass |

## Python API

```python
from peff_uniprot_fetcher import fetch_peff, fetch_peff_to_file, fasta_to_peff, fasta_to_peff_file
from pefftacular import write_peff

# Fetch and write in one call
fetch_peff_to_file("human.peff", query="organism_id:9606 AND reviewed:true")

# Or get the data back
header, entries = fetch_peff(accessions=["P12345", "Q99999"])
write_peff(header, entries, "output.peff")

# From a local FASTA file
fasta_to_peff_file("input.fasta", "output.peff")
header, entries = fasta_to_peff("input.fasta")
```

Annotation behaviour is controlled via an `AnnotationConfig` dataclass, or individual keyword arguments:

```python
from peff_uniprot_fetcher import AnnotationConfig, fetch_peff_to_file

# Using the config dataclass
cfg = AnnotationConfig(include_glycosylation=True, only_known_mass=True)
fetch_peff_to_file("human.peff", query="organism_id:9606 AND reviewed:true", cfg=cfg)

# Or pass flags directly as keyword arguments
fetch_peff_to_file("human.peff", query="organism_id:9606 AND reviewed:true", include_glycosylation=True)
```

| Parameter | Default | Effect |
|---|---|---|
| `include_variants` | `True` | Include sequence variants |
| `include_modifications` | `True` | Include PTMs (`ModResPsi`, `ModResUnimod`, `ModRes`) |
| `include_processed` | `True` | Include processed forms |
| `include_glycosylation` | `False` | Include glycosylation sites (resolved via PTM ontologies when possible) |
| `include_lipidation` | `False` | Include lipidation sites (resolved via PTM ontologies when possible) |
| `include_crosslinks` | `False` | Include cross-links (`ModRes`) |
| `only_known_mass` | `False` | Only include modifications with a known monoisotopic mass |

## PEFF annotations

The following UniProt GFF feature types are mapped to PEFF annotations:

| UniProt feature | PEFF key |
|---|---|
| Natural variant, Mutagenesis, Sequence conflict | `VariantSimple` / `VariantComplex` |
| Alternative sequence (isoform) | `VariantComplex` |
| Modified residue (PSI-MOD cross-ref) | `ModResPsi` |
| Modified residue (UniMod cross-ref) | `ModResUnimod` |
| Modified residue (UniProt PTM match, no PSI-MOD/UniMod cross-ref) | `ModRes` |
| Glycosylation (PTM match) | `ModResPsi` / `ModResUnimod` / `ModRes` |
| Lipidation (PTM match) | `ModResPsi` / `ModResUnimod` / `ModRes` |
| Cross-link | `ModRes` |
| Signal peptide | `Processed` (`PEFF:0001001`) |
| Transit peptide | `Processed` (`PEFF:0001002`) |
| Propeptide | `Processed` (`PEFF:0001003`) |
| Chain (mature protein) | `Processed` (`PEFF:0001004`) |
| Peptide | `Processed` (`PEFF:0001005`) |

A modified residue with both a PSI-MOD and UniMod cross-reference appears in both `ModResPsi` and `ModResUnimod` simultaneously. Modifications that cannot be resolved to a known PTM entry are silently skipped.

### PTM name resolution

`ModResPsi` and `ModResUnimod` entries use the **canonical ontology name** from [psimodpy](https://pypi.org/project/psimodpy/) and [unimodpy](https://pypi.org/project/unimodpy/) rather than the UniProt ptmlist name. For example, a phosphoserine site is written as `2|MOD:00046|O-phospho-L-serine` instead of `2|MOD:00046|Phosphoserine`.

PTM entries that have PSI-MOD or UniMod cross-references but lack a formula or mass in UniProt's ptmlist (loaded via [uniprotptmpy](https://pypi.org/project/uniprotptmpy/)) are automatically enriched with masses and formulas from the PSI-MOD / UniMod databases at load time.

## Human proteome script

`scripts/human_proteome_peff.py` generates a PEFF file for the reviewed human proteome and prints per-feature-type modification statistics (PSI-MOD / UniMod / both / custom / none counts, has-mass counts, top N modification names).

```bash
uv run python scripts/human_proteome_peff.py [OUTPUT] [--query QUERY] \
    [--include-glycosylation] [--include-lipidation] [--include-crosslinks] \
    [--no-variants] [--no-modifications] [--no-processed] \
    [--only-known-mass] [--top-n N]
```

`OUTPUT` defaults to `human_proteome.peff`. The `--include-*` flags opt in to feature types that are off by default; `--no-*` flags turn off features that are on by default.

## Just recipes

```bash
just download-ecoli        # download raw E. coli K-12 FASTA + GFF to data/ecoli/
just fetch-ecoli           # generate PEFF for E. coli K-12
just fasta-to-peff-ecoli   # convert downloaded E. coli FASTA to PEFF
```

## Development

```bash
just lint      # ruff check
just format    # ruff format
just check     # lint + type check + test
just test      # pytest
```

## License

[MIT](https://github.com/tacular-omics/peff_uniprot_fetcher/blob/main/LICENSE)
