"""
STANDARD_GAS_COMPOSITIONS.py

Standard gaseous mixture compositions are stored in this file in the form of
dictionaries.

The keys are the constituent species' names,
same as in the 'gas_specific_heat_coeffs_7.yaml' file (from where majority of the thermodynamic data is sourced).

The values are the mole fraction of the respective constituent in the mixture
"""

from typing import Dict

AIR: Dict[str, float] = {"N2": 0.78084, "O2": 0.209476, "Ar": 0.00934, "CO2": 0.000314}
