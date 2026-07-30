#!/bin/bash
#SBATCH --job-name=mmseqs_decontam
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --mail-user=idan.slurm@gmail.com
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/mmseqs_decontam_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/mmseqs_decontam_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

set -euo pipefail

echo "=== Starting mmseqs_decontam job ==="
date

############################################
# 0. ENVIRONMENT
############################################

. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh

conda activate TD2_env

############################################
# 0.1 TAXONOMIC ROOT TO REMOVE
############################################
# Taxonomic roots to remove (comma-separated NCBI taxids)
# Examples:
#   Bacteria: 2, Archaea: 2157, Fungi: 4751, Viridiplantae: 33090, Viruses: 10239
REMOVE_ROOTS="2,2157,4751,33090,10239"

############################################
# 1. INPUTS & PATHS
############################################

INPUT_DIR=$(readlink -f "$1")
if [[ -z "${INPUT_DIR:-}" ]]; then
    echo "No input directory provided! Exiting."
    exit 1
fi

FINAL_NAME=$(basename "$INPUT_DIR") 

shopt -s nullglob
REF_FASTA_FILE=("$INPUT_DIR"/cd-hit/*_CD-HIT.fasta)
if [[ ${#REF_FASTA_FILE[@]} -eq 0 ]]; then
    echo "CRITICAL ERROR: Could not find the clustered CD-HIT FASTA file in $INPUT_DIR/cd-hit"
    exit 1
fi
FINAL_REF="${REF_FASTA_FILE[0]}"

MMSEQ_DB_DIR="/sci/labs/ariel.chipman/idansh/mmseq_db"
TARGET_DB="${MMSEQ_DB_DIR}/uniref90"

TMP_DIR="${INPUT_DIR}/mmseqs_tmp"
OUT_DIR="${INPUT_DIR}/mmseqs_decontam"
FINAL_DIR="${INPUT_DIR}/decontamination"

mkdir -p "$TMP_DIR" "$OUT_DIR" "$FINAL_DIR"

TAX_TSV_WITH_IDS="${FINAL_DIR}/taxAssignments.with_taxid.tsv"
NAME_TAXID_MAP="${FINAL_DIR}/name_taxid_map.tsv"
QUERY_DB="${OUT_DIR}/queryDB"
TAX_DB="${OUT_DIR}/taxDB"
TAX_TSV="${FINAL_DIR}/taxAssignments.tsv"

KEEP_IDS="${FINAL_DIR}/keep_ids.txt"
CONTAM_IDS="${FINAL_DIR}/contaminants.ids"
UNCLASS_IDS="${FINAL_DIR}/unclassified.ids"

CLEAN_FASTA="${INPUT_DIR}/${FINAL_NAME}_decontam.fasta"
CONTAM_FASTA="${FINAL_DIR}/${FINAL_NAME}_contaminants.fasta"
UNCLASS_FASTA="${FINAL_DIR}/${FINAL_NAME}_unclassified.fasta"

KEPT_REPORT="${FINAL_DIR}/${FINAL_NAME}_report_kept.tsv"
TOSSED_REPORT="${FINAL_DIR}/${FINAL_NAME}_report_tossed.tsv"

echo "Input FASTA: $FINAL_REF"
echo "Target DB:   $TARGET_DB"
echo "Tmp dir:     $TMP_DIR"
echo "Out dir:     $OUT_DIR"
echo "Clean FASTA: $CLEAN_FASTA"

############################################
# 2. CREATE QUERY DB
############################################

if [[ -f "${QUERY_DB}.dbtype" ]]; then
    echo "--- Query DB ${QUERY_DB} exists, skipping createdb ---"
else
    echo "--- Creating MMseqs query DB from FASTA ---"
    mmseqs createdb "$FINAL_REF" "$QUERY_DB"
fi

############################################
# 3. TAXONOMY ASSIGNMENT
############################################

if [[ -f "${TAX_DB}.dbtype" ]]; then
    echo "--- Taxonomy DB ${TAX_DB} exists, skipping taxonomy ---"
else
    echo "--- Running mmseqs taxonomy ---"
    mmseqs taxonomy "$QUERY_DB" "$TARGET_DB" "$TAX_DB" "$TMP_DIR" \
        --threads 12 \
        --max-seqs 50 \
        --split-memory-limit 140G
fi

############################################
# 4. EXPORT TAXONOMY ASSIGNMENTS TO TSV
############################################
# This is the same call that worked for you before.
# If this step ever segfaults again, you can skip recreating it when it already exists.

if [[ -s "$TAX_TSV" ]]; then
    echo "--- Taxonomy TSV exists at $TAX_TSV, skipping createtsv ---"
else
    echo "--- Exporting taxonomy assignments to TSV ---"
    mmseqs createtsv "$QUERY_DB" "$TARGET_DB" "$TAX_DB" "$TAX_TSV"
fi

############################################
# 4.1 ADD NUMERIC TAXIDS FROM NAMES (tax_env)
############################################

if [[ -s "$TAX_TSV_WITH_IDS" ]]; then
    echo "--- Augmented TSV with taxids exists at $TAX_TSV_WITH_IDS, skipping add_taxids_from_names ---"
else
    echo "--- Adding numeric taxids based on names using tax_env ---"

    NAME_COL=3  # adjust after inspecting TAX_TSV

    conda deactivate
    conda activate tax_env

    /sci/labs/ariel.chipman/idansh/scripts/transcriptome_assembly/add_taxids_from_names.py \
        "$TAX_TSV" "$NAME_COL" "$TAX_TSV_WITH_IDS" "$NAME_TAXID_MAP"

    conda deactivate
    conda activate TD2_env
fi

############################################
# 5. BUILD KEEP LIST BY LINEAGE (REMOVE_ROOTS) - EXCLUDING UNCLASSIFIED
############################################

if [[ -s "$KEEP_IDS" ]]; then
    echo "--- Keep list exists at $KEEP_IDS, skipping build ---"
else
    echo "--- Building keep_ids.txt from augmented TSV (removing lineages under $REMOVE_ROOTS) ---"
    conda deactivate
    conda activate tax_env

    TAXID_COL=-1  # last column (numeric taxid)

    python - "$TAX_TSV_WITH_IDS" "$TAXID_COL" "$KEEP_IDS" "$REMOVE_ROOTS" << 'EOF'
import sys
from ete3 import NCBITaxa

tsv = sys.argv[1]
taxid_col = int(sys.argv[2])  # -1 = last column
keep_out = sys.argv[3]
remove_roots_str = sys.argv[4]

# Roots to remove, e.g. Bacteria=2, Fungi=4751, etc.
REMOVE_ROOTS = {int(t) for t in remove_roots_str.split(",") if t.strip().isdigit()}

# Taxids that should be considered "unclassified" - do NOT add to keep
UNCLASSIFIED_TAXIDS = {
    1,       # root
    12908,   # unclassified sequences
    32644,   # unidentified
    131567,  # cellular organisms (too vague)
}

ncbi = NCBITaxa()

keep = set()
to_resolve = set()

# First pass: collect all taxids we need lineages for
with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        col = taxid_col if taxid_col >= 0 else (len(parts) - 1)
        taxid_str = parts[col].strip()
        if taxid_str.isdigit():
            to_resolve.add(int(taxid_str))

# Precompute lineages
lineages = {}
for t in to_resolve:
    try:
        lineages[t] = set(ncbi.get_lineage(t))
    except Exception:
        lineages[t] = set()

# Second pass: decide which queries to keep (ONLY properly classified non-contaminants)
with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        col = taxid_col if taxid_col >= 0 else (len(parts) - 1)
        qid = parts[0]
        taxid_str = parts[col].strip()
        
        # Skip sequences without numeric taxids (unclassified)
        if not taxid_str.isdigit():
            continue
        
        taxid = int(taxid_str)
        
        # Skip sequences with unclassified taxids
        if taxid in UNCLASSIFIED_TAXIDS:
            continue
            
        lineage = lineages.get(taxid, set())
        # If any of the REMOVE_ROOTS is in the lineage, this is a contaminant
        if any(root in lineage for root in REMOVE_ROOTS):
            continue
        # Otherwise, keep (properly classified, non-contaminant)
        keep.add(qid)

with open(keep_out, "w") as out:
    for q in sorted(keep):
        out.write(q + "\n")
        
print(f"Added {len(keep)} sequences to keep list (excluding contaminants and unclassified)")
EOF

    conda deactivate
    conda activate TD2_env
fi

N_KEEP=$(wc -l < "$KEEP_IDS" || echo 0)
echo "Number of contigs to keep (properly classified, excluding roots $REMOVE_ROOTS and unclassified): $N_KEEP"

############################################
# 5.2 BUILD CONTAMINANT + UNCLASSIFIED ID LISTS
############################################

TAXID_COL=-1  # last column (numeric taxid)

if [[ -s "$CONTAM_IDS" && -s "$UNCLASS_IDS" ]]; then
    echo "--- Contaminant and unclassified ID lists already exist, skipping build ---"
else
    echo "--- Classifying remaining contigs as contaminants or unclassified ---"
    conda deactivate
    conda activate tax_env

    python - "$TAX_TSV_WITH_IDS" "$TAXID_COL" "$KEEP_IDS" "$CONTAM_IDS" "$UNCLASS_IDS" << 'EOF'
import sys
from ete3 import NCBITaxa

tsv = sys.argv[1]
taxid_col = int(sys.argv[2])  # -1 = last column
keep_file = sys.argv[3]
contam_out = sys.argv[4]
unclass_out = sys.argv[5]

BACTERIA = 2
ARCHAEA = 2157
FUNGI = 4751
PLANTS = 33090      # Viridiplantae
VIRUSES = 10239

CONTAM_ROOTS = {BACTERIA, ARCHAEA, FUNGI, PLANTS, VIRUSES}

# Taxids that should be considered "unclassified"
UNCLASSIFIED_TAXIDS = {
    1,       # root
    12908,   # unclassified sequences
    32644,   # unidentified
    131567,  # cellular organisms (too vague)
}

ncbi = NCBITaxa()

keep = set(x.strip() for x in open(keep_file) if x.strip())
contam = set()
unclass = set()

to_resolve = set()

with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        col = taxid_col if taxid_col >= 0 else (len(parts) - 1)
        taxid_str = parts[col].strip()
        if taxid_str.isdigit():
            to_resolve.add(int(taxid_str))

lineages = {}
for t in to_resolve:
    try:
        lineages[t] = set(ncbi.get_lineage(t))
    except Exception:
        lineages[t] = set()

with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        qid = parts[0]
        if qid in keep:
            continue
        col = taxid_col if taxid_col >= 0 else (len(parts) - 1)
        taxid_str = parts[col].strip()
        
        # Check for invalid or unclassified taxids
        if not taxid_str.isdigit():
            unclass.add(qid)
            continue
            
        taxid = int(taxid_str)
        
        # Check if taxid is in the unclassified set
        if taxid in UNCLASSIFIED_TAXIDS:
            unclass.add(qid)
            continue
            
        lineage = lineages.get(taxid, set())
        if any(root in lineage for root in CONTAM_ROOTS):
            contam.add(qid)
        else:
            unclass.add(qid)

with open(contam_out, "w") as out:
    for q in sorted(contam):
        out.write(q + "\n")

with open(unclass_out, "w") as out:
    for q in sorted(unclass):
        out.write(q + "\n")
EOF

    conda deactivate
    conda activate TD2_env
fi

N_CONTAM=$(wc -l < "$CONTAM_IDS" || echo 0)
N_UNCLASS=$(wc -l < "$UNCLASS_IDS" || echo 0)
echo "Contaminant sequences: $N_CONTAM"
echo "Unclassified sequences: $N_UNCLASS"

############################################
# 5.1 BUILD REPORTS FOR KEPT AND TOSSED CONTIGS
############################################

echo "--- Building kept/tossed reports ---"
python - "$FINAL_REF" "$TAX_TSV" "$KEEP_IDS" "$KEPT_REPORT" "$TOSSED_REPORT" << 'EOF'
import sys

fasta = sys.argv[1]
tsv = sys.argv[2]
keep_file = sys.argv[3]
kept_report = sys.argv[4]
tossed_report = sys.argv[5]

TAXID_COL = 1
RANK_COL = 2
NAME_COL = 3

keep = set(x.strip() for x in open(keep_file) if x.strip())

tax_info = {}
with open(tsv) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= max(TAXID_COL, NAME_COL):
            continue
        qid = parts[0]
        taxid = parts[TAXID_COL]
        rank = parts[RANK_COL] if len(parts) > RANK_COL else ""
        name = parts[NAME_COL] if len(parts) > NAME_COL else ""
        tax_info[qid] = (taxid, rank, name)

lengths = {}
with open(fasta) as f:
    current_id = None
    current_len = 0
    for line in f:
        if line.startswith(">"):
            if current_id is not None:
                lengths[current_id] = current_len
            current_id = line[1:].strip().split()[0]
            current_len = 0
        else:
            current_len += len(line.strip())
    if current_id is not None:
        lengths[current_id] = current_len

def write_report(out_path, ids):
    ids = [i for i in ids if i in lengths]
    total = len(ids)
    if total == 0:
        with open(out_path, "w") as out:
            out.write("id\tlength\ttaxid\trank\tname\n")
            out.write("# No entries\n")
        return
    sizes = [lengths[i] for i in ids]
    avg_len = sum(sizes) / total
    min_len = min(sizes)
    max_len = max(sizes)
    with open(out_path, "w") as out:
        out.write("id\tlength\ttaxid\trank\tname\n")
        for qid in ids:
            taxid, rank, name = tax_info.get(qid, ("NA", "NA", "NA"))
            out.write(f"{qid}\t{lengths[qid]}\t{taxid}\t{rank}\t{name}\n")
        out.write(f"# total\t{total}\n")
        out.write(f"# avg_length\t{avg_len:.2f}\n")
        out.write(f"# min_length\t{min_len}\n")
        out.write(f"# max_length\t{max_len}\n")

all_ids = set(lengths.keys())
kept_ids = sorted(keep & all_ids)
tossed_ids = sorted(all_ids - keep)

write_report(kept_report, kept_ids)
write_report(tossed_report, tossed_ids)
EOF

echo "Kept report:    $KEPT_REPORT"
echo "Tossed report:  $TOSSED_REPORT"

############################################
# 6. WRITE THREE FASTAS (+ unclassified in clean)
############################################

echo "--- Writing FASTAs for kept/contaminant/unclassified ---"
python - "$FINAL_REF" "$KEEP_IDS" "$CONTAM_IDS" "$UNCLASS_IDS" \
        "$CLEAN_FASTA" "$CONTAM_FASTA" "$UNCLASS_FASTA" << 'EOF'
import sys

fasta = sys.argv[1]
keep_file = sys.argv[2]
contam_file = sys.argv[3]
unclass_file = sys.argv[4]
fa_kept = sys.argv[5]
fa_contam = sys.argv[6]
fa_unclass = sys.argv[7]

def load_ids(path):
    return set(x.strip() for x in open(path) if x.strip())

keep = load_ids(keep_file)
contam = load_ids(contam_file)
unclass = load_ids(unclass_file)

# Open all output files
out_keep = open(fa_kept, "w")
out_contam = open(fa_contam, "w")
out_unclass = open(fa_unclass, "w")

with open(fasta) as f:
    current_id = None
    current_category = None
    
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            current_id = name
            
            if name in keep:
                current_category = "keep"
                out_keep.write(line)
            elif name in contam:
                current_category = "contam"
                out_contam.write(line)
            elif name in unclass:
                current_category = "unclass"
                # Write to BOTH clean and unclassified FASTAs
                out_keep.write(line)
                out_unclass.write(line)
            else:
                current_category = None
        else:
            # Write sequence lines to appropriate file(s)
            if current_category == "keep":
                out_keep.write(line)
            elif current_category == "contam":
                out_contam.write(line)
            elif current_category == "unclass":
                # Write to BOTH clean and unclassified FASTAs
                out_keep.write(line)
                out_unclass.write(line)

out_keep.close()
out_contam.close()
out_unclass.close()

print(f"Wrote {len(keep)} kept sequences to {fa_kept}")
print(f"Wrote {len(contam)} contaminant sequences to {fa_contam}")
print(f"Wrote {len(unclass)} unclassified sequences to {fa_unclass}")
print(f"Clean FASTA contains {len(keep) + len(unclass)} total sequences (kept + unclassified)")
EOF

rm -rf "$OUT_DIR" "$TMP_DIR"

echo "Kept FASTA (includes unclassified):    $CLEAN_FASTA"
echo "Contaminants FASTA:                     $CONTAM_FASTA"
echo "Unclassified FASTA (separate tracking): $UNCLASS_FASTA"
