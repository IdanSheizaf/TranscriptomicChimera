#!/usr/bin/env python
import sys
from ete3 import NCBITaxa

if len(sys.argv) != 4:
    sys.stderr.write("Usage: build_keep_ids.py taxAssignments.tsv TAXID_COL keep_ids.txt\n")
    sys.exit(1)

tsv = sys.argv[1]
taxid_col = int(sys.argv[2])  # 0-based index of taxid column
keep_out = sys.argv[3]

ARTHROPODA_TAXID = 6656

ncbi = NCBITaxa()

keep = set()
to_resolve = set()

# First pass: collect all taxids
with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= taxid_col:
            continue
        taxid_str = parts[taxid_col].strip()
        if not taxid_str.isdigit():
            continue
        to_resolve.add(int(taxid_str))

# Precompute lineages for all taxids seen
lineages = {}
for t in to_resolve:
    try:
        lineages[t] = set(ncbi.get_lineage(t))
    except Exception:
        lineages[t] = set()

# Second pass: decide which queries to keep
with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= taxid_col:
            continue
        qid = parts[0]
        taxid_str = parts[taxid_col].strip()
        if not taxid_str.isdigit():
            continue
        taxid = int(taxid_str)
        if ARTHROPODA_TAXID in lineages.get(taxid, set()):
            keep.add(qid)

with open(keep_out, "w") as out:
    for q in sorted(keep):
        out.write(q + "\n")
