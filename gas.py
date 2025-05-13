import math
from typing import Dict, Set, List, Tuple
from periodictable.formulas import formula
from cache.CONSTANTS import UNIVERSAL_GAS_CONSTANT_SI as R_SI


class GasMixture:
    def __init__(
        self,
        name: str,
        mole_fraction_composition: Dict[str, float],
        thermo_data_yaml: Dict[str, Dict[str, Dict]],
        least_admissible_mole_frac_error: float = 1e-4,
    ):
        """Create a GasMixture object to model the behaviour of mixture of ideal gas species. Stores the parsed data from NASA_THERMO_DATA_YAML,
        which includes valid temperature ranges and their respective 7 coefficients

        Args:
            name (str): Name of the gas mixture
            mole_fraction_composition (Dict[str, float]): Dictionary of gas mixture composition where the key is the constituent species and
            the value is the mole fraction of species in the mixture
            least_admissible_mole_frac_error (float): Error threshold to validate composition
            thermo_data_yaml (Dict[str, Dict[str, Dict]]): User side pre-loaded dictionary of thermal data
        """
        # name of the mixture
        self.__name: str = name

        # set of all the species present in the mixture
        self.__species: Set[str] = set(mole_fraction_composition.keys())

        # Constituent species mole fraction
        self.__mole_fraction_composition: Dict[str, float] = mole_fraction_composition

        # Ensure all the mole fractions add up to ~ 1
        self.__validate_molar_composition(least_admissible_mole_frac_error)

        # Availabe temperature ranges and 7 coefficients stored for each species
        self.__species_data: Dict[str, Dict] = {}

        # Load the data from the NASA THERMO DATA YAML file
        self.__load_constituent_species_data(thermo_data_yaml)

    @property
    def mole_fraction_composition(self):
        return self.__mole_fraction_composition

    @property
    def species(self):
        return self.__species

    @property
    def name(self):
        return self.__name

    @property
    def species_data(self):
        return self.__species_data

    def __load_constituent_species_data(self, yaml_data: Dict[str, Dict[str, Dict]]):
        """Stores the temperature range and
        coefficients in the self._species_data dictionary for each constituent species

        Args:
            yaml_data (Dict[str, Dict[str, Dict]]): User side pre-loaded dictionary of thermal data

        Raises:
            ValueError: In case of any species whose data is not availabe in yaml_data file
        """
        species_map = {
            str(species["name"]).lower(): species["thermo"]
            for species in yaml_data["species"]
        }

        for constituent_species in self.__species:
            key = constituent_species.lower()
            if key not in species_map:
                raise ValueError(
                    f"{self.__name}: Invalid species: '{constituent_species}'"
                )
            self.__species_data[constituent_species] = species_map[key]
            self.__species_data[constituent_species]["molar-mass"] = (
                (formula(constituent_species).molecular_mass) * 1e-3 * 6.022e23
            )

            # Remove irrelevant fields safely
            self.__species_data[constituent_species].pop("model", None)
            self.__species_data[constituent_species].pop("note", None)

    def __validate_molar_composition(self, least_admissible_error: float):
        """Ensures all the mole fractions of the constituent species add up to 1 within permitted range of error

        Args:
            least_admissible_error (float): Minimum threshold error above which molar composition is invalid

        Raises:
            ValueError: If all the mole fractions do not add up to 1 within the least_admissible_error range.
        """
        total = sum(self.__mole_fraction_composition.values())
        if not abs(total - 1.0) < least_admissible_error:
            raise ValueError(
                f"{self.__name}: Total mole fraction must sum to 1.0. Got: {total}"
            )

    def __str__(self):
        lines = [f"MOLE FRACTION COMPOSITION FOR '{self.__name}':\n"]

        for k, v in self.__mole_fraction_composition.items():
            lines.append(f"{k}: {v * 100:.6f} %")

        lines.append(f"\nSPECIES DATA FOR '{self.__name}':\n")

        for k, v in self.__species_data.items():
            lines.append(f"'{k}'\n")
            lines.append(f"Temperature Ranges: {v['temperature-ranges']}")
            if len(v["temperature-ranges"]) == 3:
                lines.append(f"Lower Range Coefficients: {v['data'][0]}")
                lines.append(f"Higher Range Coefficients: {v['data'][1]}\n")
            else:
                lines.append(f"Coefficients: {v['data'][0]}\n")

        return "\n".join(lines)


class GasState:
    def __init__(self, gas_mixture: GasMixture, temperature_K: float = 298.0):
        """Stores all the derived thermodynamic data calculated from the stored coefficients in a Gas Mixture

        Args:
            gas_mixture (GasMixture): GasMixture object that stores all the relevant coefficients for thermodynamic calculations
            temperature_K (float): Gas temperature, which is the independent variable for specific heats, enthalpy and entropy functions.
            Defaults to 298. (Standard Air Temperature).
        """
        # Gas mixture to be evaluated
        self.__gas_mixture: GasMixture = gas_mixture

        # Gas temperature
        self.__temperature_K: float = temperature_K

        # Isobaric molar specific heat
        self.__cp_SI: float = None

        # Isochoric molar specific heat
        self.__cv_SI: float = None

        # Molar enthalpy
        self.__enthalpy_SI: float = None

        # Molar entropy
        self.__entropy_SI: float = None

        # Ratio of specific heats
        self.__gamma: float = None

        # Molar mass of the mixture in kg
        self.__mol_mass_kg: float = None

        # Calculates all the thermodynamic properties for a given temperature
        self.__calc_thermodynamic_properties()

    def __thermo_nasa_polynomials(
        self, coefficients: List[float]
    ) -> Tuple[float, float, float]:
        """NASA polynomial function which uses the 7 coefficients to calculate
        the isobaric specific heat, molar enthalpy and molar entropy of each species

        Args:
            coefficients (List[float]): Nasa coefficients

        Returns:
            Tuple[float, float, float]: (isobaric specific heat/R, enthalpy /RT, entropy/R)
        """
        t = self.__temperature_K
        a1, a2, a3, a4, a5, a6, a7 = coefficients
        cp_R = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4
        h_RT = (
            a1
            + (a2 * t / 2)
            + (a3 * t**2 / 3)
            + (a4 * t**3 / 4)
            + (a5 * t**4 / 5)
            + (a6 / t)
        )

        s_R = (
            a1 * math.log(t)
            + a2 * t
            + (a3 / 2) * t**2
            + (a4 / 3) * t**3
            + (a5 / 4) * t**4
            + a7
        )
        return cp_R, h_RT, s_R

    def __calc_thermodynamic_properties(self):
        """Calculates isobaric specific heat, molar enthalpy and molar entropy for the gas mixture

        Raises:
            ValueError: if the gas temperature is out of the range for which data is available
        """
        cp_total = 0.0
        enthalpy_total = 0.0
        entropy_total = 0.0
        molar_mass_total = 0.0
        gas_mixture_species = self.__gas_mixture.species

        for species in gas_mixture_species:
            species_data = self.__gas_mixture.species_data[species]
            temp_range = species_data["temperature-ranges"]
            coeff_array = species_data["data"]
            species_mol_mass_kg = species_data["molar-mass"]
            temp_range_len = len(temp_range)

            if temp_range_len == 3:
                t_lo, t_mid, t_hi = temp_range

                if t_lo <= self.__temperature_K < t_mid:
                    coeffs = coeff_array[0]

                elif t_mid <= self.__temperature_K <= t_hi:
                    coeffs = coeff_array[1]

                else:
                    raise ValueError(
                        f"[{species}] T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            else:
                t_lo, t_hi = temp_range

                if t_lo <= self.__temperature_K <= t_hi:
                    coeffs = coeff_array[0]
                else:
                    raise ValueError(
                        f"[{species}] T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            cp_R, h_RT, s_R = self.__thermo_nasa_polynomials(coeffs)
            mole_frac = self.__gas_mixture.mole_fraction_composition[species]
            cp = cp_R * mole_frac
            enthalpy = h_RT * mole_frac
            entropy = s_R * mole_frac
            molar_mass = mole_frac * species_mol_mass_kg

            cp_total += cp
            enthalpy_total += enthalpy
            entropy_total += entropy
            molar_mass_total += molar_mass

        self.__cp_SI = cp_total * R_SI
        self.__enthalpy_SI = enthalpy_total * R_SI * self.__temperature_K
        self.__entropy_SI = entropy_total * R_SI
        self.__cv_SI = self.__cp_SI - R_SI
        self.__gamma = self.__cp_SI / self.__cv_SI
        self.__mol_mass_kg = molar_mass_total

    def updateTemperature(self, temperature_K: float):
        """Update the gas temperature and recalculate the thermodynamic properties

        Args:
            temperature_K (float): Updated gas temperature
        """
        self.__temperature_K = temperature_K
        self.__calc_thermodynamic_properties()

    def updateGasMixture(self, gas_mixture: GasMixture):
        self.__gas_mixture = gas_mixture
        self.__calc_thermodynamic_properties()

    @property
    def cp(self):
        return self.__cp_SI / self.molar_mass

    @property
    def cv(self):
        return self.__cv_SI / self.molar_mass

    @property
    def enthalpy(self):
        return self.__enthalpy_SI / self.molar_mass

    @property
    def entropy(self):
        return self.__entropy_SI / self.molar_mass

    @property
    def gamma(self):
        return self.__gamma

    @property
    def molar_mass(self):
        return self.__mol_mass_kg

    @property
    def temperature(self):
        return self.__temperature_K

    def __str__(self):
        lines = [self.__gas_mixture.__str__()]
        lines.append(f"\nTHERMODYNAMIC DATA FOR '{self.__gas_mixture.name}':\n")
        lines.append(f"Molar Mass (kg): {self.molar_mass}")
        lines.append(f"Temperature (K): {self.temperature}")
        lines.append(f"Isobaric specific heat  (J / kg K): {self.cp}")
        lines.append(f"Isochoric specific heat (J / kg K): {self.cv}")
        lines.append(f"Specific Enthalpy (J / kg): {self.enthalpy}")
        lines.append(f"Specific Entropy (J / kg K): {self.entropy}")
        lines.append(f"Ratio of specific heats: {self.gamma}")

        return "\n".join(lines)
