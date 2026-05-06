"""
01_setup_environment.py
=======================
Environment & path bootstrap for the AFDB-truncated-dimer pipeline.

Mounts Google Drive (Colab), creates the project base directory, and
prepares the working directory so subsequent steps (data download,
mapping, profiling, figures) can use plain relative file names.

Run order: 1 of 14
"""

import os
from pathlib import Path

import pandas as pd
import polars as pl  # Essential for 31M rows

from google.colab import drive

# 1. Mount Google Drive
drive.mount('/content/drive')

# 2. Define and create the project base path
BASE_PATH = Path('/content/drive/MyDrive/MIT/All Project Files/AFDB dimer')
BASE_PATH.mkdir(parents=True, exist_ok=True)

# 3. Make the base path the working directory.
# In a Colab cell this is normally done with the "%cd" magic:
#     %cd "/content/drive/MyDrive/MIT/All Project Files/AFDB dimer"
# When run as a plain Python script, fall back to os.chdir():
os.chdir(BASE_PATH)

print(f"Current Working Directory: {os.getcwd()}")
