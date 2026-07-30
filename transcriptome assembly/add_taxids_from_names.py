#!/usr/bin/env python
import sys
from ete3 import NCBITaxa

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: add_taxids_from_names.py taxAssignments.tsv NAME_COL OUT_TSV OUT_MAP\n"
    )
    sys.exit(1)

tsv_in = sys.argv[1]
name_col = int(sys.argv[2])  # 0-based column index for taxon name/lineage
tsv_out = sys.argv[3]
map_out = sys.argv[4]        # optional: name -> taxid map for inspection

ncbi = NCBITaxa()

names = set()
rows = []

with open(tsv_in) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            rows.append((line, None))
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= name_col:
            rows.append((line, None))
            continue
        name = parts[name_col].strip()
        rows.append((line, name))
        if name:
            names.add(name)

# Resolve names -> taxids
name_to_taxid = {}
if names:
    # get_name_translator expects a list of names and returns a dict
    trans = ncbi.get_name_translator(list(names))
    for nm in names:
        # Some names may not resolve; skip them
        taxids = trans.get(nm, [])
        name_to_taxid[nm] = taxids[0] if taxids else None

# Write augmented TSV with extra taxid column at the end
with open(tsv_out, "w") as out:
    for line, name in rows:
        if not line.strip() or line.startswith("#") or name is None:
            out.write(line)
            continue
        taxid = name_to_taxid.get(name)
        taxid_str = str(taxid) if taxid is not None else "NA"
        line = line.rstrip("\n")
        out.write(f"{line}\t{taxid_str}\n")

# Optional: write mapping table for debugging
with open(map_out, "w") as out:
    out.write("name\ttaxid\n")
    for nm in sorted(names):
        taxid = name_to_taxid.get(nm)
        out.write(f"{nm}\t{taxid if taxid is not None else 'NA'}\n")
