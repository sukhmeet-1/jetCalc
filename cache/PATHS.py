"""
PATHS.py

Folder and file paths frequently used in the project are cached here.
"""

import os


PROJECT_ROOT_PATH: str = os.getcwd()
DATA_PATH: str = os.path.join(PROJECT_ROOT_PATH, "data")
DATA_NASA_PATH: str = os.path.join(DATA_PATH, "nasa")

# Thermodynamic data file used for calculation of gas properties
NASA_THERMAL_DATA_PATH: str = os.path.join(DATA_NASA_PATH, "nasa_thermal_data.yaml")
