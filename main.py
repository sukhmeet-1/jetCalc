import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import AIR
from mixture import GasMixture

# TODO parse polynomial coefficients for cp and cv calculation of mixture of gaseous species
# TODO write a function for calculating cp and cv for a given input temperature


def main():
    with open(COEFF_7_NASA_CP_DATA_PATH, "r") as f:
        data = yaml.safe_load(f)

    print(
        f"Thermodynamic Data available for {len(data['species'])} gaseous species, stored at {COEFF_7_NASA_CP_DATA_PATH}"
    )

    # An example of how mixtures of gaseous species can be constructed

    air = GasMixture("Standard Air", AIR)
    air.mixture_info()


if __name__ == "__main__":
    main()
