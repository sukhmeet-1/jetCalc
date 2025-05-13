import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from mixture import GasMixture

# TODO write a function for calculating cp and cv for a given input temperature


def main():
    with open(COEFF_7_NASA_CP_DATA_PATH, "r") as f:
        data = yaml.safe_load(f)

    print(
        f"Thermodynamic Data available for {len(data['species'])} gaseous species, stored at {COEFF_7_NASA_CP_DATA_PATH}"
    )

    # An example of how mixtures of gaseous species can be constructed

    air = GasMixture(name="Standard Air", mole_fraction_composition=STANDARD_AIR_COMPOSITION)
    air.mixture_info()


if __name__ == "__main__":
    main()
