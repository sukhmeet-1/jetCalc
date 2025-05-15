from typing import Dict, Set
from periodictable.formulas import formula
from cache.CONSTANTS import AVOGADRO_NUM as AVOG_N


class GasMixture:
    def __init__(
        self,
        name: str,
        mole_fraction_composition: Dict[str, float],
        thermo_data_yaml: Dict[str, Dict[str, Dict]],
        admissible_mole_frac_error: float = 1e-4,
    ):
        """Initializes a mixture of gaseous species by reading thermodynamic data from a yaml file.

        Args:
            name (str): Name of the mixture
            mole_fraction_composition (Dict[str, float]): Dictionary of constituent gaseous species and their respective mole fractions
            thermo_data_yaml (Dict[str, Dict[str, Dict]]): Thermal data loaded from a yaml file already parsed on the user side
            admissible_mole_frac_error (float, optional): Maximum error in mole fraction composition (making sure all species mole fractions add up to 1). Defaults to 1e-4.
        """
        # name of the mixture
        self.__name: str = name

        # set of all the species present in the mixture
        self.__species: Set[str] = set(mole_fraction_composition.keys())

        # Constituent species mole fraction
        self.__mole_fraction_composition: Dict[str, float] = mole_fraction_composition

        # Ensure all the mole fractions add up to ~ 1
        self.__validate_molar_composition(admissible_mole_frac_error)

        # Availabe temperature ranges and 7 coefficients stored for each species
        self.__species_data: Dict[str, Dict] = {}

        # Load the parsed species data
        self.__load_constituent_species_data(thermo_data_yaml)

    def __load_constituent_species_data(self, yaml_data: Dict[str, Dict[str, Dict]]):
        """Loads the thermal data for each of the constituent species of the gas mixture. Data includes valid temperature range and the respective 7 NASA thermal coefficients.

        Args:
            yaml_data (Dict[str, Dict[str, Dict]]): Thermal data from the yaml file parsed on the user side

        Raises:
            ValueError: If the gas mixture has a constituent species which is either unavailable or is tagged by a different name in the yaml file
        """
        # For checking if the all constituent species are available in the file
        species_map: Dict[str, Dict] = {
            str(species["name"]).lower(): species["thermo"]
            for species in yaml_data["species"]
        }

        for constituent_species in self.__species:
            key: str = constituent_species.lower()

            if key not in species_map:
                raise ValueError(
                    f"{self.__name}: Invalid species: '{constituent_species}'"
                )

            self.__species_data[constituent_species] = species_map[key]
            self.__species_data[constituent_species]["molar-mass"] = (
                (formula(constituent_species).molecular_mass) * 1e-3 * AVOG_N
            )  # molecular mass property gives the mass of 1 molecule. Multiplied with 0.001 and Avogadro's number to get mass of 1 mole in kg

            # Remove irrelevant fields safely
            self.__species_data[constituent_species].pop("model", None)
            self.__species_data[constituent_species].pop("note", None)

    def __validate_molar_composition(self, admissible_error: float):
        """Validates if all the mass fractions of the constituent species add up to 1

        Args:
            admissible_error (float): Maximum allowed error in the difference between the actual sum of mole fractions and 1

        Raises:
            ValueError: If the sum of all the mole fracions of the constituent species do not add up to 1 within the admissible error range
        """
        total: float = sum(self.__mole_fraction_composition.values())
        if not abs(total - 1.0) < admissible_error:
            raise ValueError(
                f"{self.__name}: Total mole fraction must sum to 1.0. Got: {total}"
            )

    @property
    def mole_fraction_composition(self) -> Dict[str, float]:
        return self.__mole_fraction_composition

    @property
    def species(self) -> Set[str]:
        return self.__species

    @property
    def name(self) -> str:
        return self.__name

    @property
    def species_data(self) -> Dict[str, Dict]:
        return self.__species_data

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
