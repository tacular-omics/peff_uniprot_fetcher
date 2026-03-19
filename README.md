# peff_uniprot_fetcher

Generate [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format) files from the UniProt REST API. Fetches protein sequences and GFF annotations (variants, PTMs, processed forms) and writes them as annotated PEFF using [pefftacular](https://github.com/pgarrett-scripps/pefftacular).

## Installation

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

| Flag | Effect |
|---|---|
| `--no-variants` | Exclude sequence variants (`VariantSimple`, `VariantComplex`) |
| `--no-modifications` | Exclude PTMs (`ModResPsi`, `ModRes`) |
| `--no-processed` | Exclude processed forms (`Signal peptide`, `Chain`, etc.) |

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

## PEFF annotations

The following UniProt GFF feature types are mapped to PEFF annotations:

| UniProt feature | PEFF key |
|---|---|
| Natural variant, Mutagenesis, Sequence conflict | `VariantSimple` / `VariantComplex` |
| Alternative sequence (isoform) | `VariantComplex` |
| Modified residue | `ModResPsi` (with PSI-MOD accession) or `ModRes` |
| Glycosylation, Lipidation, Cross-link | `ModRes` |
| Signal peptide | `Processed` (`PEFF:0001001`) |
| Transit peptide | `Processed` (`PEFF:0001002`) |
| Propeptide | `Processed` (`PEFF:0001003`) |
| Chain (mature protein) | `Processed` (`PEFF:0001004`) |
| Peptide | `Processed` (`PEFF:0001005`) |

PTM names are mapped to PSI-MOD accessions using the [UniProt PTM list](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt).

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
