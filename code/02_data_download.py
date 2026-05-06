"""
02_data_download.py
===================
Probe the EBI AFDB FTP collaboration tree, then download:
  (a) the AFDB multimer manifest (model_entity_metadata_mapping.csv, ~4.58 GB)
  (b) the full AFDB sequences.fasta (~118 GB)

Both downloads use aria2 with 16 parallel connections to keep wall time
manageable.  T

In Colab the shell magics ("%cd", "!apt-get", "!aria2c") run the embedded
shell commands.

Run order: 2 of 14
"""

from ftplib import FTP
import os
from pathlib import Path

import polars as pl


# ---------------------------------------------------------------------------
# 1. Probe the AFDB / NVDA FTP directory to confirm where the manifest lives
# ---------------------------------------------------------------------------
def scan_nvda_collab():
    """List the AFDB / NVIDIA collaboration directory on the EBI FTP."""
    try:
        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()
        # The 2026 quaternary expansion is hosted here
        path = '/pub/databases/alphafold/collaborations/nvda/'
        ftp.cwd(path)

        items = []
        ftp.retrlines('LIST', items.append)
        print(f"--- Contents of {path} ---")
        for item in items:
            print(item)

        ftp.quit()
    except Exception as e:
        print(f"Error: {e}")


# ---------------------------------------------------------------------------
# 2. Download the manifest (Colab shell magics)
# ---------------------------------------------------------------------------
# %cd "/content/drive/MyDrive/MIT/All Project Files/AFDB dimer"
# !apt-get install -y aria2 > /dev/null
# !aria2c -x 16 -s 16 \
#     "ftp://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/model_entity_metadata_mapping.csv"


# ---------------------------------------------------------------------------
# 3. Ingest the manifest in Polars (uses RAM heavily) and filter for human
#    homodimers based on taxId and ipTM availability.
# ---------------------------------------------------------------------------
def ingest_manifest():
    BASE_PATH = Path('/content/drive/MyDrive/MIT/All Project Files/AFDB dimer')
    os.chdir(BASE_PATH)

    print("Ingesting 31M complex records...")
    df = pl.read_csv("model_entity_metadata_mapping.csv", ignore_errors=True)

    # 3a. Human only (TaxID 9606)
    human_complexes = df.filter(pl.col("taxId") == 9606)

    # 3b. Restrict to entries with ipTM (i.e. are actually complexes)
    homodimers = human_complexes.filter(pl.col("ipTM").is_not_null())
    print(f"Human Complexes identified: {len(homodimers)}")
    return homodimers


# ---------------------------------------------------------------------------
# 4. Build a fast Human-only FASTA lookup (only run if you need to subset
#    sequences.fasta -- the 31M-row workflow is faster with the full lookup
#    in 05_state_transition_profiling.py).
# ---------------------------------------------------------------------------
def build_human_sequence_lookup(fasta_path):
    print("Building Human Sequence Lookup (RAM-intensive)...")
    seq_map = {}
    current_acc = None
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # >tr|A0A024RBG1|A0A024RBG1_HUMAN ...
                current_acc = (line.split('|')[1] if '|' in line
                               else line[1:].split()[0])
            elif current_acc:
                seq_map[current_acc] = seq_map.get(current_acc, "") + line.strip()
    return seq_map


# ---------------------------------------------------------------------------
# 5. Download the full 118 GB sequences.fasta (only if missing)
# ---------------------------------------------------------------------------
# %cd "/content/drive/MyDrive/MIT/All Project Files/AFDB dimer"
# !apt-get install -y aria2 > /dev/null
# !aria2c -x 16 -s 16 -c \
#     "ftp://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta"


if __name__ == "__main__":
    scan_nvda_collab()
    # ingest_manifest()
    # seq_map = build_human_sequence_lookup("sequences.fasta")
