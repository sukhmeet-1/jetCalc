import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from mixture import GasMixture, GasMixtureThermoProperties


def main():
    with open(COEFF_7_NASA_CP_DATA_PATH, "r") as f:
        data = yaml.safe_load(f)

    print(
        f"Thermodynamic Data available for {len(data['species'])} gaseous species, stored at {COEFF_7_NASA_CP_DATA_PATH}\n"
    )

    # An example of how mixtures of gaseous species can be constructed

    air = GasMixture(
        name="Standard Air", mole_fraction_composition=STANDARD_AIR_COMPOSITION
    )
    airThermo = GasMixtureThermoProperties(gas_mixture=air, temperature_K=350)
    print(airThermo)


if __name__ == "__main__":
    main()
