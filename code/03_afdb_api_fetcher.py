"""
03_afdb_api_fetcher.py
======================
Per-protein AlphaFold Database (AFDB) fetcher.

For a given UniProt canonical accession, this script:
  1. Pulls the canonical sequence from the UniProt REST API.
  2. Verifies stoichiometry-2 (homodimer) entries via the AFDB API and
     downloads PDB + JSON metadata for the canonical.
  3. Walks UniProt 'ALTERNATIVE PRODUCTS' isoforms, keeps only those
     that are a proper subset of the canonical sequence (natural
     truncations), and downloads their AFDB models too.

This is the "single-protein" companion to the manifest-driven pipeline
in 02_data_download.py / 04_homodimer_mapping.py.  It is most useful for
building targeted small-cohort datasets around a curated list of leads.

Run order: 3 of 14
"""

import json
import time
from pathlib import Path

import requests

from google.colab import drive


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
drive.mount('/content/drive')
BASE_PATH = Path('/content/drive/MyDrive/MIT/All Project Files/AFDB dimer')
BASE_PATH.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def is_natural_truncation(canonical_seq, isoform_seq):
    """Verify the isoform is a strict, shorter substring of the canonical."""
    if not canonical_seq or not isoform_seq:
        return False
    return isoform_seq in canonical_seq and len(isoform_seq) < len(canonical_seq)


def fetch_and_save_afdb(uniprot_id, category_path):
    """Fetch a homodimer PDB + JSON from the AFDB API. Skip if already on disk."""
    pdb_path = category_path / f"{uniprot_id}.pdb"
    json_path = category_path / f"{uniprot_id}.json"

    if pdb_path.exists() and json_path.exists():
        return True  # already downloaded

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        time.sleep(0.1)  # polite to the AFDB API
        response = requests.get(api_url, timeout=10)
        if response.status_code != 200:
            return False

        data = response.json()
        # Keep only stoichiometry == 2 (homodimer)
        dimer_entries = [e for e in data if e.get('stoichiometry') in [[2], 2]]
        if not dimer_entries:
            return False

        entry = dimer_entries[0]  # most recent
        pdb_url = entry['pdbUrl']

        with open(pdb_path, 'wb') as f:
            f.write(requests.get(pdb_url).content)
        with open(json_path, 'w') as f:
            json.dump(entry, f)
        return True
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return False


def process_protein_cluster(canonical_id):
    """Map one canonical to its truncated isoforms and fetch all dimer models."""
    cluster_dir = BASE_PATH / canonical_id
    if cluster_dir.exists():
        pass  # placeholder for incremental logic

    try:
        u_url = f"https://rest.uniprot.org/uniprotkb/{canonical_id}?format=json"
        u_res = requests.get(u_url)
        if u_res.status_code != 200:
            return
        u_data = u_res.json()
        canonical_seq = u_data.get('sequence', {}).get('value', '')

        # 1. Canonical dimer
        canon_path = cluster_dir / 'canonical'
        canon_path.mkdir(parents=True, exist_ok=True)
        if not fetch_and_save_afdb(canonical_id, canon_path):
            return

        # 2. Walk isoforms
        isoforms = []
        for comment in u_data.get('comments', []):
            if comment.get('commentType') == 'ALTERNATIVE PRODUCTS':
                for iso in comment.get('isoforms', []):
                    for acc in iso.get('isoformIds', []):
                        if acc != canonical_id:
                            isoforms.append(acc)

        for iso_id in isoforms:
            i_url = f"https://rest.uniprot.org/uniprotkb/{iso_id}?format=json"
            i_data = requests.get(i_url).json()
            iso_seq = i_data.get('sequence', {}).get('value', '')

            if is_natural_truncation(canonical_seq, iso_seq):
                iso_path = cluster_dir / f"isoform_{iso_id}"
                iso_path.mkdir(parents=True, exist_ok=True)
                fetch_and_save_afdb(iso_id, iso_path)

    except Exception as e:
        print(f"Cluster {canonical_id} failed: {e}")


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Replace with your real list (or load from a CSV column 'Entry')
    target_list = ["P12345", "Q9Y2I1"]
    for p_id in target_list:
        process_protein_cluster(p_id)
        print(f"Finished Cluster: {p_id}")
