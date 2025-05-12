import yaml
from paths import *

# TODO parse polynomial coefficients for cp and cv calculation of mixture of gaseous species


def main(species_name):
    with open(COEFF_7_NASA_CP_DATA_PATH, "r") as f:
        data = yaml.safe_load(f)

    print(f"cp / R polynomial coefficients available for the following species:")
    for species in data["species"]:
        print(species["name"])


if __name__ == "__main__":
    main("CO2")
