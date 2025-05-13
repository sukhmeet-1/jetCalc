import math
from typing import Dict, Set, List
import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH as NASA_THERMO_DATA_YAML
from cache.CONSTANTS import UNIVERSAL_GAS_CONSTANT_SI


class GasMixture:
    def __init__(
        self,
        name: str,
        mole_fraction_composition: Dict[str, float],
    ):
        self._name: str = name

        self._species: Set[str] = set(mole_fraction_composition.keys())
        self._mole_fraction_composition: Dict[str, float] = mole_fraction_composition
        self._validate_molar_composition()

        self._species_data: Dict[str, Dict] = {}
        self._load_constituent_species_data()

        self._species_temperature_ranges: Dict[str, List[float]] = {}
        self._species_coefficients_array: Dict[str, List[List[float]]] = {}
        self._load_constituent_temperature_ranges_and_coefficients()

    def _load_constituent_species_data(self):
        with open(NASA_THERMO_DATA_YAML, "r") as f:
            data = yaml.safe_load(f)

        available_species_names = {
            str(species["name"]).lower() for species in data["species"]
        }
        invalid_species = [
            s for s in self._species if s.lower() not in available_species_names
        ]
        if invalid_species:
            raise ValueError(f"Invalid species: {invalid_species}")

        constituent_species_data_dict = {}
        for constituent_species in self._species:
            for species in data["species"]:
                if str(species["name"]).lower() == constituent_species.lower():
                    constituent_species_data_dict[constituent_species] = species[
                        "thermo"
                    ]
                    break
        self._species_data = constituent_species_data_dict

    def _load_constituent_temperature_ranges_and_coefficients(self):
        for species in self._species:
            self._species_temperature_ranges[species] = self._species_data[species][
                "temperature-ranges"
            ]
            self._species_coefficients_array[species] = self._species_data[species][
                "data"
            ]

    def _validate_molar_composition(self, least_admissible_error=1e-4):
        total = sum(self._mole_fraction_composition.values())
        if not abs(total - 1.0) < least_admissible_error:
            raise ValueError(f"Total mole fraction must sum to 1.0. Got: {total}")

    def info(self):
        print(f"Mole Fraction Composition of {self._name} is:")

        for k, v in self._mole_fraction_composition.items():
            print(f"{k}: {v*100}%")
        print()

        for k, v in self._species_data.items():
            print(k)
            print(v)
            print()


class GasMixtureThermoProperties:
    def __init__(self, gas_mixture: GasMixture, temperature_K=298):
        self.gas_mixture = gas_mixture
        self.temperature_K = temperature_K
        self.cp_SI = 0
        self.cv_SI = 0
        self.enthalpy_SI = 0
        self.entropy_SI = 0
        self.gamma = 0

        self._calc_thermo(temperature_K)

    def _thermo_polynomials(self, temperature_K, coefficients: List[float]):
        a1, a2, a3, a4, a5, a6, a7 = coefficients
        cp_R = (
            a1
            + a2 * temperature_K
            + a3 * math.pow(temperature_K, 2)
            + a4 * math.pow(temperature_K, 3)
            + a5 * math.pow(temperature_K, 4)
        )
        h_RT = (
            a1
            + (a2 / 2) * temperature_K
            + (a3 / 3) * math.pow(temperature_K, 2)
            + (a4 / 4) * math.pow(temperature_K, 3)
            + (a5 / 5) * math.pow(temperature_K, 4)
            + (a6 / temperature_K)
        )

        s_R = (
            a1 * math.log(temperature_K)
            + a2 * temperature_K
            + (a3 / 2) * math.pow(temperature_K, 2)
            + (a4 / 3) * math.pow(temperature_K, 3)
            + (a5 / 4) * math.pow(temperature_K, 4)
            + a7
        )
        return cp_R, h_RT, s_R

    def _calc_thermo(self, temperature_K):
        cp_total = 0.0
        enthalpy_total = 0.0
        entropy_total = 0.0

        for species in self.gas_mixture._species:
            temp_range_len = len(self.gas_mixture._species_temperature_ranges[species])

            if temp_range_len == 3:
                t_lo, t_mid, t_hi = self.gas_mixture._species_temperature_ranges[
                    species
                ]
                if t_lo <= temperature_K < t_mid:
                    coeffs = self.gas_mixture._species_coefficients_array[species][0]

                elif t_mid <= temperature_K <= t_hi:
                    coeffs = self.gas_mixture._species_coefficients_array[species][1]

                else:
                    raise ValueError(
                        f"Temperature {temperature_K}K out of range for species '{species}'"
                    )

            else:
                t_lo, t_hi = self.gas_mixture._species_temperature_ranges[species]
                if t_lo <= temperature_K <= t_hi:
                    coeffs = self.gas_mixture._species_coefficients_array[species][0]
                else:
                    raise ValueError(
                        f"Temperature {temperature_K}K out of range for species '{species}'"
                    )

            cp_R, h_RT, s_R = self._thermo_polynomials(temperature_K, coeffs)
            mole_frac = self.gas_mixture._mole_fraction_composition[species]
            cp = cp_R * mole_frac
            enthalpy = h_RT * mole_frac
            entropy = s_R * mole_frac

            cp_total += cp
            enthalpy_total += enthalpy
            entropy_total += entropy
        self.cp_SI = cp_total * UNIVERSAL_GAS_CONSTANT_SI
        self.enthalpy_SI = enthalpy_total * UNIVERSAL_GAS_CONSTANT_SI * temperature_K
        self.entropy_SI = entropy_total * UNIVERSAL_GAS_CONSTANT_SI
        self.cv_SI = self.cp_SI - UNIVERSAL_GAS_CONSTANT_SI
        self.gamma = self.cp_SI / self.cv_SI

    def update_temperature(self, temperature_K: float):
        self._calc_thermo(temperature_K)

    def info(self):
        self.gas_mixture.info()
        print("THERMODYNAMIC DATA\n")
        print(f"Temperature (K): {self.temperature_K}")
        print(f"Isobaric molar specific heat  (SI): {self.cp_SI}")
        print(f"Isochoric molar specific heat (SI): {self.cv_SI}")
        print(f"Molar Enthalpy (SI): {self.enthalpy_SI}")
        print(f"Molar Entropy (SI): {self.entropy_SI}")
        print(f"Ratio of specific heats: {self.gamma}")


if __name__ == "__main__":
    pass
