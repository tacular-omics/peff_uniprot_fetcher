# Copilot instructions

The canonical guide for AI agents in this repo is [`CLAUDE.md`](../CLAUDE.md) at the repo
root: commands, module layout, public API, conventions and gotchas. Read it first.

Key rules:

1. Use `just` recipes (`just check` = ruff + ty + pytest), falling back to `uv run`.
   Never use pip or call `python` / `pytest` outside `uv run`.
2. Tests must stay offline. Never add a test that calls the UniProt REST API; build
   FASTA/GFF strings inline instead.
3. `_client.py` must keep working under Pyodide (the GitHub Pages app in `docs/`):
   keep `httpx` import-guarded and the fetch path synchronous.
4. Do not bump the version, tag, or publish. Only the tacular-omics overseer releases.
5. Python >= 3.12 typing (`X | None`, builtin generics), frozen dataclasses for config,
   short docstrings, `logging` not `print` in `src/`.
