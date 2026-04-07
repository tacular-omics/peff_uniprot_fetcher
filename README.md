# peff_uniprot_fetcher

Generate [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format) files from the UniProt REST API. Fetches protein sequences and GFF annotations (variants, PTMs, processed forms) and writes them as annotated PEFF using [pefftacular](https://github.com/pgarrett-scripps/pefftacular).

> **Try it in the browser — no install required:** <https://tacular-omics.github.io/peff_uniprot_fetcher/>

## Web app

A static, zero-backend web app is hosted on GitHub Pages. Type an NCBI taxonomy ID (e.g. `83333` for *E. coli* K-12, `9606` for human), click **Generate PEFF**, and the browser fetches the UniProt data, builds an annotated PEFF file, and hands you a download link — all without touching a server.

It runs the same `peff_uniprot_fetcher` Python package you'd use from the CLI, compiled to WebAssembly via [Pyodide](https://pyodide.org) inside a Web Worker. UniProt is called directly from the browser, so the page is fully client-side.

- **Open it:** <https://tacular-omics.github.io/peff_uniprot_fetcher/>
- **When to use the CLI instead:** large proteomes (human, mouse, plants) pull hundreds of MB of GFF and take several minutes in-tab — use the CLI or Python API below for those. The web app is best for small/medium organisms and quick one-offs.
- **Source:** the static bundle lives in [`docs/`](docs/) and loads the project wheel via `micropip`; see `docs/worker.js` for the boot sequence.

## Installation

Only needed if you want the CLI or Python API — skip this if you're using the web app above.

```bash
# From source
git clone https://github.com/pgarrett-scripps/peff_uniprot_fetcher
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
| Modified residue (UniProt PTM match) | `ModRes` |
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

`ModResPsi` and `ModResUnimod` entries use the **canonical ontology name** from [tacular](https://github.com/pgarrett-scripps/tacular) (PSI-MOD or UniMod) rather than the UniProt ptmlist name. For example, a phosphoserine site is written as `2|MOD:00046|O-phospho-L-serine` instead of `2|MOD:00046|Phosphoserine`.

PTM entries that have PSI-MOD or UniMod cross-references but lack a formula or mass in UniProt's ptmlist are automatically enriched with masses and formulas from tacular at load time.

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
