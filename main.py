import yaml
from cache.DIRECTORIES import NASA_THERMAL_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from thermodynamics.phase import GasMixture
from thermodynamics.state import GasState

with open(NASA_THERMAL_DATA_PATH, "r") as f:
    NASA_THERMAL_DATA = yaml.safe_load(f)


def main():
    air = GasMixture(
        name="Standard Air",
        mole_fraction_composition=STANDARD_AIR_COMPOSITION,
        thermo_data_yaml=NASA_THERMAL_DATA,
    )
    state = GasState(
        name="Random", gas_mixture=air, pressure_Pa=101325, density_SI=1.225
    )
    print(state)
    state.update_state( pressure_Pa=202325, density_SI=1.225)
    print(state)


if __name__ == "__main__":
    main()
