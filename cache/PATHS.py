"""
PATHS.py

Folder and file paths frequently used in the project are cached here.
"""

import os


PROJECT_ROOT_PATH: str = os.getcwd()
DATA_PATH: str = os.path.join(PROJECT_ROOT_PATH, "data")
DATA_NASA_PATH: str = os.path.join(DATA_PATH, "nasa")

# Thermodynamic data file used for calculation of gas properties
COEFF_7_NASA_CP_DATA_PATH: str = os.path.join(
    DATA_NASA_PATH, "gas_specific_heat_coeffs_7.yaml"
)
