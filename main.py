import yaml
from cache.PATHS import COEFF_7_NASA_CP_DATA_PATH as NASA_THERMAL_DATA
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from mixture import GasMixture, GasMixtureThermoProperties


def main():
    with open(NASA_THERMAL_DATA, "r") as f:
        nasa_gas_thermal_data = yaml.safe_load(f)

    air = GasMixture(
        name="Standard Air",
        mole_fraction_composition=STANDARD_AIR_COMPOSITION,
        thermo_data_yaml=nasa_gas_thermal_data,
    )
    airThermoProperties = GasMixtureThermoProperties(gas_mixture=air, temperature_K=350)
    print(airThermoProperties)


if __name__ == "__main__":
    main()
