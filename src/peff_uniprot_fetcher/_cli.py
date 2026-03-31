"""CLI entry points for peff_uniprot_fetcher."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from peff_uniprot_fetcher import AnnotationConfig, fasta_to_peff_file, fetch_peff_to_file
from peff_uniprot_fetcher._client import fetch_entry, stream_search


def _configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)


def _annotation_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--no-variants", action="store_true", help="Exclude sequence variants.")
    parser.add_argument("--no-modifications", action="store_true", help="Exclude PTM/modification annotations.")
    parser.add_argument(
        "--no-processed",
        action="store_true",
        help="Exclude processed form annotations (signal peptides, chains, etc.).",
    )
    parser.add_argument(
        "--glycosylation", action="store_true", help="Include glycosylation annotations (default: off)."
    )
    parser.add_argument("--lipidation", action="store_true", help="Include lipidation annotations (default: off).")
    parser.add_argument("--crosslinks", action="store_true", help="Include cross-link annotations (default: off).")
    parser.add_argument(
        "--only-known-mass",
        action="store_true",
        help="Only include modifications with a known monoisotopic mass.",
    )


def _annotation_config(args: argparse.Namespace) -> AnnotationConfig:
    return AnnotationConfig(
        include_variants=not args.no_variants,
        include_modifications=not args.no_modifications,
        include_processed=not args.no_processed,
        include_glycosylation=args.glycosylation,
        include_lipidation=args.lipidation,
        include_crosslinks=args.crosslinks,
        only_known_mass=args.only_known_mass,
    )


def fasta_to_peff_cli() -> None:
    """Entry point: fasta-to-peff"""
    parser = argparse.ArgumentParser(
        prog="fasta-to-peff",
        description="Convert a local UniProt FASTA file to PEFF by fetching GFF annotations from UniProt.",
    )
    parser.add_argument("input", metavar="INPUT", help="Path to input UniProt FASTA file.")
    parser.add_argument("output", metavar="OUTPUT", help="Path for the output PEFF file.")
    _annotation_flags(parser)
    args = parser.parse_args()
    _configure_logging()

    fasta_to_peff_file(fasta=args.input, output=args.output, cfg=_annotation_config(args))
    logging.info("Done. Written to %s", args.output)


def fetch_peff_cli() -> None:
    """Entry point: fetch-peff"""
    parser = argparse.ArgumentParser(
        prog="fetch-peff",
        description="Fetch proteins from UniProt and write a PEFF file.",
    )
    parser.add_argument("output", metavar="OUTPUT", help="Path for the output PEFF file.")

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--organism-id",
        "-x",
        metavar="TAXID",
        help="NCBI taxonomy ID (e.g. 9606 for human). Fetches all reviewed (Swiss-Prot) entries.",
    )
    source.add_argument(
        "--query",
        "-q",
        metavar="QUERY",
        help="Raw UniProt search query (e.g. 'organism_id:9606 AND reviewed:true').",
    )
    source.add_argument(
        "--accessions",
        "-a",
        nargs="+",
        metavar="ACC",
        help="One or more UniProt accessions (e.g. P12345 Q99999).",
    )

    parser.add_argument(
        "--unreviewed",
        action="store_true",
        help="When using --organism-id, include unreviewed (TrEMBL) entries (default: reviewed only).",
    )
    _annotation_flags(parser)
    args = parser.parse_args()
    _configure_logging()

    cfg = _annotation_config(args)

    def _call(*, query: str | None = None, accessions: list[str] | None = None) -> None:
        fetch_peff_to_file(output=args.output, query=query, accessions=accessions, cfg=cfg)

    if args.organism_id:
        reviewed_clause = "" if args.unreviewed else " AND reviewed:true"
        _call(query=f"organism_id:{args.organism_id}{reviewed_clause}")
    elif args.query:
        _call(query=args.query)
    else:
        _call(accessions=args.accessions)

    logging.info("Done. Written to %s", args.output)


def download_uniprot_cli() -> None:
    """Entry point: download-uniprot"""
    parser = argparse.ArgumentParser(
        prog="download-uniprot",
        description="Download raw FASTA and/or GFF files from UniProt for inspection.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default=".",
        metavar="DIR",
        help="Directory to write files into (default: current directory).",
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--organism-id", "-x", metavar="TAXID", help="NCBI taxonomy ID (e.g. 9606 for human).")
    source.add_argument("--query", "-q", metavar="QUERY", help="Raw UniProt search query.")
    source.add_argument("--accession", "-a", metavar="ACC", help="Single UniProt accession.")

    parser.add_argument(
        "--unreviewed",
        action="store_true",
        help="When using --organism-id, include unreviewed entries.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["fasta", "gff"],
        choices=["fasta", "gff"],
        metavar="FMT",
        help="Formats to download: fasta, gff (default: both).",
    )
    args = parser.parse_args()
    _configure_logging()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.organism_id:
        reviewed_clause = "" if args.unreviewed else " AND reviewed:true"
        query = f"organism_id:{args.organism_id}{reviewed_clause}"
        stem = f"uniprot_organism_{args.organism_id}"
        for fmt in args.formats:
            logging.info("Downloading %s...", fmt.upper())
            text = stream_search(query, fmt=fmt)
            path = out_dir / f"{stem}.{fmt}"
            path.write_text(text)
            logging.info("Written %d bytes to %s", len(text), path)

    elif args.query:
        stem = "uniprot_query"
        for fmt in args.formats:
            logging.info("Downloading %s...", fmt.upper())
            text = stream_search(args.query, fmt=fmt)
            path = out_dir / f"{stem}.{fmt}"
            path.write_text(text)
            logging.info("Written %d bytes to %s", len(text), path)

    else:
        acc = args.accession
        for fmt in args.formats:
            logging.info("Downloading %s for %s...", fmt.upper(), acc)
            text = fetch_entry(acc, fmt=fmt)
            path = out_dir / f"{acc}.{fmt}"
            path.write_text(text)
            logging.info("Written %d bytes to %s", len(text), path)
