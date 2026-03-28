default: lint format check test

# Install dependencies
install:
    uv sync

# Run linting checks
lint:
    uv run ruff check src

# Format code
format:
	uv run ruff check --select I --fix src
	uv run ruff format src

# Run ty type checker
ty:
    uv run ty check src

# Run type checking
check:
    just lint
    just ty
    just test

# Run tests
test:
    uv run pytest tests

# Download raw E. coli K-12 FASTA and GFF from UniProt for inspection (tax ID 83333)
download-ecoli:
    uv run download-uniprot --organism-id 83333 --output-dir data/ecoli

# Generate PEFF for E. coli K-12 (reviewed Swiss-Prot entries only)
fetch-ecoli:
    uv run fetch-peff data/ecoli/ecoli.peff --organism-id 83333

# Convert a locally downloaded E. coli FASTA to PEFF
fasta-to-peff-ecoli fasta="data/ecoli/uniprot_organism_83333.fasta":
    uv run fasta-to-peff {{fasta}} data/ecoli/ecoli_from_fasta.peff

# Remove build artifacts
clean:
    rm -rf dist

# Build the package and check that the .obo data file is included
build:
    uv build
