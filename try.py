import sys
from pefftacular import read_peff

path = sys.argv[1] if len(sys.argv) > 1 else "human.peff"

print(f"Reading {path}...")
header, entries = read_peff(path)

print(f"PEFF version: {header.peff_version}")
for db in header.databases:
    print(f"  {db.db_name}: {db.number_of_entries} entries")

print(f"\nTotal entries read: {len(entries)}")

annotated = [e for e in entries if e.variant_simple or e.mod_res_psi or e.processed]
print(f"Entries with annotations: {len(annotated)}")

for entry in annotated[:3]:
    print(f"\n  {entry.prefix}:{entry.db_unique_id} {entry.pname}")
    if entry.variant_simple:
        print(f"    VariantSimple:  {entry.variant_simple[:3]}{'...' if len(entry.variant_simple) > 3 else ''}")
    if entry.mod_res_psi:
        print(f"    ModResPsi:      {entry.mod_res_psi[:3]}{'...' if len(entry.mod_res_psi) > 3 else ''}")
    if entry.processed:
        print(f"    Processed:      {entry.processed}")
