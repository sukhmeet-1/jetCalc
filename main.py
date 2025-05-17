import yaml
from cache.PATHS import NASA_THERMAL_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from thermodynamics.phase import GasMixture
from thermodynamics.state import GasState
from thermodynamics.process import (
    GasProcess,
    GasProcessConstraint,
    GasProcessType,
    generate_process_constraints,
)

with open(NASA_THERMAL_DATA_PATH, "r") as f:
    NASA_THERMAL_DATA = yaml.safe_load(f)


def main():
    air = GasMixture(
        name="Standard Air",
        mole_fraction_composition=STANDARD_AIR_COMPOSITION,
        thermo_data_yaml=NASA_THERMAL_DATA,
    )
    state_1 = GasState(
        gas_mixture=air, mass_kg=3, pressure_Pa=101325, temperature_K=340
    )
    const: GasProcessConstraint = generate_process_constraints(
        volume=state_1.volume * 3, mass=state_1.mass * 0.9
    )
    process = GasProcess(GasProcessType.ISOTHERMAL, state_1, const)
    print(process.initial_state)
    print(process.final_ideal_state)
    print(
        f"Ideal work done in the {process.process_type} process: {process.ideal_work_done}"
    )


if __name__ == "__main__":
    main()
