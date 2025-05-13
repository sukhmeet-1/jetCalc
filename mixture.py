from typing import Dict, Set


class GasMixture:
    def __init__(self, name: str, mole_fraction_dict: Dict[str, float]):
        self._name: str = name
        self._species: Set[str] = set(mole_fraction_dict.keys())
        self._mole_fraction_composition: Dict[str, float] = mole_fraction_dict
        self.cp: float = -1.0
        self.cv: float = -1.0
        self.gamma: float = -1.0

    def mixture_info(self):
        print(
            f"Mole Fraction Composition of {self._name} is {self._mole_fraction_composition}"
        )


if __name__ == "__main__":
    pass
