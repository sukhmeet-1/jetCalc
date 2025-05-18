from typing import Dict, Set, List, Any
from periodictable.formulas import formula
from cache.CONSTANTS import AVOGADRO_NUM as AVOG_N


class GasMixture:
    def __init__(
        self,
        name: str = None,
        mole_fraction_composition: Dict[str, float] = None,
        thermo_data_yaml: Dict[str, Any] = None,
        admissible_mole_frac_error: float = 1e-6,
    ):
        """Initializes a mixture of gaseous species by using external thermodynamic data.

        Args:
            `name` (str): Name of the mixture
            `mole_fraction_composition` (Dict[str, float]): Dictionary of constituent gaseous species and their respective mole fractions
            `thermo_data_yaml` (Dict[str, Any]): External thermal data of gaseous species
            `admissible_mole_frac_error` (float, optional): Maximum error in mole fraction composition (making sure all species mole fractions add up to 1). Defaults to 1e-4.
        Raises:
            ValueError: If the gas mixture is not named
            ValueError: If the mole fraction composition of the mixture is not provided or is empty
            ValueError: If the constituent species do not have their names capitalized
            ValueError: If the mole fraction of any of the constituent species is zero or negative
            ValueError: If the thermal data parsed from a yaml file is not provided
            ValueError: If the sum of all the mole fractions of the constituent species do not add up to 1 within the admissible error range
            ValueError: If the parsed yaml data has an invalid structure - 'species', 'thermo', 'temperature-ranges' or 'data' keys are missing
            ValueError: If the any one of the constituent species is not available in the parsed yaml data
        """

        if name is None:
            raise ValueError("Name must be provided for the gas mixture")

        if mole_fraction_composition is None or mole_fraction_composition == {}:
            raise ValueError(
                "Mole fraction composition of the gas mixture must be provided"
            )

        for species_names in mole_fraction_composition.keys():
            if species_names.lower() == species_names:
                raise ValueError("The species names must be all capitalized")

        for species, mole_fraction in mole_fraction_composition.items():
            if mole_fraction <= 0:
                raise ValueError(
                    f"Mole fraction of the constituent species, {species} is {mole_fraction}. It should be positive and non-zero"
                )

        if thermo_data_yaml is None:
            raise ValueError("Parsed yaml thermal data of the species must be provided")

        # name of the mixture
        self.__name: str = name

        # set of all the species present in the mixture
        self.__species: Set[str] = set(mole_fraction_composition.keys())

        # Constituent species mole fraction
        self.__mole_fraction_composition: Dict[str, float] = mole_fraction_composition

        # Ensure all the mole fractions add up to ~ 1
        self.__validate_molar_composition(admissible_mole_frac_error)

        # Availabe temperature ranges and 7 coefficients stored for each species
        self.__species_data: Dict[str, Any] = {}

        # Load the parsed species data
        self.__load_constituent_species_data(thermo_data_yaml)

    def __load_constituent_species_data(self, yaml_data: Dict[str, Any]):
        """Loads the thermal data for each of the constituent species of the gas mixture. Data includes valid temperature range and the respective 7 NASA thermal coefficients.

        Args:
            yaml_data (Dict[str, Any]): Thermal data from the yaml file parsed on the user side

        Raises:
            ValueError: If the parsed yaml data has an invalid structure - 'species', 'thermo', 'temperature-ranges' or 'data' keys are missing
            ValueError: If the any one of the constituent species is not available in the parsed yaml data
        """
        if "species" not in yaml_data:
            raise ValueError(
                "Provided yaml data has invalid structure. 'species' key missing"
            )

        if not isinstance(yaml_data["species"], list):
            raise ValueError("Expected 'species' to be a list of dictionaries")

        for species_data in yaml_data["species"]:
            if "thermo" not in species_data:
                raise ValueError(
                    f"Provided yaml data has invalid structure. 'thermo' key missing"
                )
            if "temperature-ranges" not in species_data["thermo"]:
                raise ValueError(
                    f"Provided yaml data has invalid structure. 'temperature-ranges' key missing"
                )
            if "data" not in species_data["thermo"]:
                raise ValueError(
                    f"Provided yaml data has invalid structure. 'data' key, containing coefficients missing"
                )
        # For checking if the all constituent species are available in the file
        species_map: Dict[str, Any] = {
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
        """Validates if all the mole fractions of the constituent species add up to 1

        Args:
            admissible_error (float): Maximum allowed relative error in the actual sum of mole fractions and 1

        Raises:
            ValueError: If the sum of all the mole fractions of the constituent species do not add up to 1 within the admissible error range
        """
        total: float = sum(self.__mole_fraction_composition.values())
        if not abs(1 - (1.0 / total)) < admissible_error:
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
    def species_data(
        self,
    ) -> Dict[str, Any]:
        return self.__species_data

    def __str__(self) -> str:
        lines: List[str] = [f"MOLE FRACTION COMPOSITION FOR '{self.__name}':\n"]

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
