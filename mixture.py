from typing import Dict, Set
import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH


class GasMixture:
    def __init__(self, name: str, mole_fraction_dict: Dict[str, float]):
        self._name: str = name
        self._species: Set[str] = set(mole_fraction_dict.keys())
        self._mole_fraction_composition: Dict[str, float] = mole_fraction_dict
        self._species_data: Dict[str, Dict] = self._load_constituent_species_data()
        self.cp: float = -1.0
        self.cv: float = -1.0
        self.gamma: float = -1.0

    def _load_constituent_species_data(self):
        with open(COEFF_7_NASA_CP_DATA_PATH, "r") as f:
            data = yaml.safe_load(f)
        constituent_species_data_dict = {}

        for constituent_species in self._species:
            for species in data["species"]:
                if str(species["name"]).lower()== constituent_species.lower():
                    constituent_species_data_dict[constituent_species] = species[
                        "thermo"
                    ]
                    break
        return constituent_species_data_dict

    def mixture_info(self):

        print(f"Mole Fraction Composition of {self._name} is:")

        for k, v in self._mole_fraction_composition.items():
            print(f"{k}: {v*100}%")
        print()
        for k, v in self._species_data.items():
            print(k)
            print(v)
            print()


if __name__ == "__main__":
    pass
