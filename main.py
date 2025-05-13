import yaml
from cache.PATHS import NASA_THERMAL_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from gas import GasMixture, GasState

with open(NASA_THERMAL_DATA_PATH, "r") as f:
    NASA_THERMAL_DATA = yaml.safe_load(f)


def main():

    air = GasMixture(
        name="Standard Air",
        mole_fraction_composition=STANDARD_AIR_COMPOSITION,
        thermo_data_yaml=NASA_THERMAL_DATA,
    )
    state = GasState(gas_mixture=air, temperature_K=298)
    print(state)


if __name__ == "__main__":
    main()
