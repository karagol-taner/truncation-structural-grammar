"""
01_setup_environment.py
=======================
Environment & path bootstrap for the AFDB-truncated-dimer pipeline.

Creates the project base directory and prepares the working directory so
subsequent steps (data download, mapping, profiling, figures) can use plain
relative file names.

Run order: 1 of 14
"""

import os
from pathlib import Path

import pandas as pd
import polars as pl

try:
    from google.colab import drive
except ImportError:
    drive = None

if drive is not None:
    drive.mount('/content/drive')

BASE_PATH = Path(os.environ.get("AFDB_DIMER_DIR", "."))
BASE_PATH.mkdir(parents=True, exist_ok=True)

os.chdir(BASE_PATH)

print(f"Current Working Directory: {os.getcwd()}")
