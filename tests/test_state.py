import pytest
import yaml
from cache.DIRECTORIES import NASA_THERMAL_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from thermodynamics.state import GasState, GasMixture

with open(NASA_THERMAL_DATA_PATH, "r") as f:
    NASA_THERMAL_DATA = yaml.safe_load(f)


@pytest.fixture(scope="module")
def standard_air_initialization():
    mixture = GasMixture(
        "Standard Air",
        STANDARD_AIR_COMPOSITION,
        NASA_THERMAL_DATA,
        admissible_mole_frac_error=1e-4,
    )
    return mixture


def state_type_validity(state: GasState):
    assert isinstance(state, GasState)
    assert isinstance(state.pressure, float)
    assert isinstance(state.volume, float)
    assert isinstance(state.mass, float)
    assert isinstance(state.temperature, float)
    assert isinstance(state.molar_mass, float)
    assert isinstance(state.density, float)
    assert isinstance(state.cp, float)
    assert isinstance(state.cv, float)
    assert isinstance(state.gamma, float)
    assert isinstance(state.R, float)
    assert isinstance(state.R_uni, float)
    assert state.R_uni == pytest.approx(8.314, rel=1e-6)
    assert isinstance(state.gas_mixture, GasMixture)
    assert isinstance(state.enthalpy, float)
    assert isinstance(state.entropy, float)


@pytest.mark.parametrize(
    "standard_input, standard_expected",
    [
        (
            {"mass_kg": 1, "volume_SI": 1, "temperature_K": 300},
            {
                "R": 287.05,
                "pressure": 86115,
                "cp": 1004.80163,
                "cv": 717.751402,
                "gamma": 1.39992988,
            },
        ),
        (
            {"pressure_Pa": 1, "volume_SI": 1, "temperature_K": 300},
            {
                "R": 287.05,
                "mass": 1.1612369e-5,
                "cp": 1004.80163,
                "cv": 717.751402,
                "gamma": 1.39992988,
            },
        ),
        (
            {"pressure_Pa": 1, "mass_kg": 1, "temperature_K": 300},
            {
                "R": 287.05,
                "volume": 86115,
                "cp": 1004.80163,
                "cv": 717.751402,
                "gamma": 1.39992988,
            },
        ),
        (
            {"pressure_Pa": 101325, "volume_SI": 1, "mass_kg": 1},
            {
                "R": 287.05,
                "temperature": 352.98699,
                "cp": 1008.883278,
                "cv": 721.833046,
                "gamma": 1.397668,
            },
        ),
        (
            {
                "pressure_Pa": 101325,
                "volume_SI": 1,
                "mass_kg": 1,
                "temperature_K": 352.98699,
            },
            {
                "R": 287.05,
                "cp": 1008.883278,
                "cv": 721.833046,
                "gamma": 1.397668,
            },
        ),
    ],
    ids=[
        "mass-volume-temp",
        "pressure-volume-temp",
        "pressure-mass-temp",
        "pressure-volume-mass",
        "all-four-specified",
    ],
)
def test_gas_state_initialization(
    standard_air_initialization, standard_input, standard_expected
):
    state = GasState(gas_mixture=standard_air_initialization, **standard_input)
    state_type_validity(state)
    for key, expected in standard_expected.items():
        assert getattr(state, key) == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize(
    "invalid_input",
    [
        # Missing mass
        {"pressure_Pa": 101325, "temperature_K": 352.98699},
        # Missing temperature
        {"pressure_Pa": 101325, "volume_SI": 1},
        # Missing pressure
        {"mass_kg": 1, "temperature_K": 352.98699},
        # Missing pressure
        {"volume_SI": 1, "mass_kg": 1},
        # Missing volume
        {"pressure_Pa": 101325, "mass_kg": 1},
        # Only pressure
        {"pressure_Pa": 101325},
        # Only temperature
        {"temperature_K": 352.98699},
        # Only mass
        {"mass_kg": 1},
        # Empty
        {},
        # negative PT
        {"pressure_Pa": -101325, "temperature_K": 352.98699},
        {"pressure_Pa": -101325, "temperature_K": -352.98699},
        {"pressure_Pa": -101325, "temperature_K": -352.98699},
        # negative PV
        {"pressure_Pa": -101325, "volume_SI": 1},
        {"pressure_Pa": -101325, "volume_SI": -1},
        {"pressure_Pa": 101325, "volume_SI": -1},
        # negative mT
        {"mass_kg": -1, "temperature_K": 352.98699},
        {"mass_kg": -1, "temperature_K": -352.98699},
        {"mass_kg": 1, "temperature_K": -352.98699},
        # negative Vm
        {"volume_SI": -1, "mass_kg": 1},
        {"volume_SI": -1, "mass_kg": -1},
        {"volume_SI": 1, "mass_kg": -1},
        # negative Pm
        {"pressure_Pa": -101325, "mass_kg": 1},
        {"pressure_Pa": -101325, "mass_kg": -1},
        {"pressure_Pa": 101325, "mass_kg": -1},
        # Only pressure
        {"pressure_Pa": -101325},
        # Only temperature
        {"temperature_K": -352.98699},
        # Only mass
        {"mass_kg": -1},
        # Temp higher than range
        {"pressure_Pa": 1, "volume_SI": 1, "temperature_K": 8000},
        # Temp lower than range
        {"pressure_Pa": 1, "volume_SI": 1, "temperature_K": 8},
    ],
    ids=[
        "pressure-temp",
        "pressure-volume",
        "mass-temp",
        "volume-mass",
        "pressure-mass",
        "pressure",
        "temperature",
        "mass",
        "no-input",
        # negative PT
        "negative-pressure-temp",
        "pressure-negative-temp",
        "negative-pressure-negative-temp",
        # negative PV
        "negative-pressure-volume",
        "pressure-negative-volume",
        "negative-pressure-negative-volume",
        # negative mT
        "negative-mass-temp",
        "mass-negative-temp",
        "negative-mass-negative-volume-temp",
        # negative Vm
        "negative-volume-mass",
        "volume-negative-mass",
        "negative-volume-negative-mass",
        # negative Pm
        "negative-pressure-mass",
        "pressure-negative-mass",
        "negative-pressure-negative-mass",
        "negative-pressure",
        "negative-temp",
        "negative-mass",
        "high-temp-range",
        "low-temp-range",
    ],
)
def test_invalid_state_initialization(standard_air_initialization, invalid_input):
    with pytest.raises(ValueError):
        GasState(gas_mixture=standard_air_initialization, **invalid_input)
