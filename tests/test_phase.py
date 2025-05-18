import pytest
import yaml
from cache.DIRECTORIES import NASA_THERMAL_DATA_PATH
from cache.STANDARD_GAS_COMPOSITIONS import STANDARD_AIR_COMPOSITION
from thermodynamics.phase import GasMixture


with open(NASA_THERMAL_DATA_PATH, "r") as f:
    NASA_THERMAL_DATA = yaml.safe_load(f)

INVALID_DATA_MISSING_SPECIES = {
    "not_species": [  # should be "species"
        {
            "name": "O2",
            "thermo": {
                "temperature-ranges": [200, 1000, 6000],
                "data": [[1, 1, 1, 1, 1, 1, 1]],
            },
        }
    ]
}

INVALID_DATA_MISSING_THERMO = {
    "species": [
        {
            "name": "O2",
            "not-thermo": {
                "temperature-ranges": [200, 1000, 6000],
                "data": [[1, 1, 1, 1, 1, 1, 1]],
            },
        }
    ]
}

INVALID_DATA_MISSING_TEMP_RANGES = {
    "species": [
        {
            "name": "O2",
            "thermo": {"data": [[1, 1, 1, 1, 1, 1, 1]]},
        }
    ]
}

INVALID_DATA_MISSING_COEFFS = {
    "species": [
        {
            "name": "O2",
            "thermo": {
                "temperature-ranges": [200, 1000, 6000],
            },
        }
    ]
}


INVALID_DATA_LIST = [
    INVALID_DATA_MISSING_SPECIES,
    INVALID_DATA_MISSING_THERMO,
    INVALID_DATA_MISSING_TEMP_RANGES,
    INVALID_DATA_MISSING_COEFFS,
]


def get_invalid_yaml_data():
    return INVALID_DATA_LIST


@pytest.fixture(scope="module")
def standard_air_initialization():
    mixture = GasMixture(
        "Standard Air",
        STANDARD_AIR_COMPOSITION,
        NASA_THERMAL_DATA,
        admissible_mole_frac_error=1e-4,
    )
    return mixture


def test_gas_mixture_initialization(standard_air_initialization):
    mixture = standard_air_initialization
    assert isinstance(mixture, GasMixture)
    assert isinstance(mixture.name, str)
    assert isinstance(mixture.species, set)
    assert isinstance(mixture.mole_fraction_composition, dict)
    for mole_frac in mixture.mole_fraction_composition.values():
        assert isinstance(mole_frac, float)
    assert isinstance(mixture.species_data, dict)

    assert mixture.name == "Standard Air"
    assert abs(sum(mixture.mole_fraction_composition.values()) - 1.0) < 1e-4
    assert set(mixture.species_data.keys()) == set(STANDARD_AIR_COMPOSITION.keys())

    molar_masses_kg = {
        "N2": 0.0280134,
        "O2": 0.0319988,
        "Ar": 0.039948,
        "CO2": 0.0440095,
        "CH4": 0.0160425,
    }

    for species in mixture.species:
        assert "temperature-ranges" in mixture.species_data[species]
        assert "data" in mixture.species_data[species]
        assert "molar-mass" in mixture.species_data[species]

        assert isinstance(mixture.species_data[species]["temperature-ranges"], list)
        assert len(mixture.species_data[species]["temperature-ranges"]) in [2, 3]
        for temp in mixture.species_data[species]["temperature-ranges"]:
            assert isinstance(temp, float)

        assert isinstance(mixture.species_data[species]["data"], list)
        assert len(mixture.species_data[species]["data"]) in [1, 2]
        for data in mixture.species_data[species]["data"]:
            assert len(data) == 7
        for array in mixture.species_data[species]["data"]:
            for coeff in array:
                assert isinstance(coeff, float)

        assert mixture.species_data[species]["molar-mass"] == pytest.approx(
            molar_masses_kg[species], rel=1e-4
        )

        assert "model" not in mixture.species_data[species]
        assert "note" not in mixture.species_data[species]


def test_invalid_composition():
    invalid_composition = {"N2": 0.7, "O2": 0.1}
    with pytest.raises(ValueError):
        GasMixture("Invalid composition", invalid_composition, NASA_THERMAL_DATA)


def test_invalid_species():
    invalid_species = {"N2": 0.7, "X3": 0.3}
    with pytest.raises(ValueError):
        GasMixture("Invalid species", invalid_species, NASA_THERMAL_DATA)


def test_zero_or_negative_mole_fractions():
    invalid_compositions = [{"O2": 0.0, "N2": 1.0}, {"O2": -0.1, "N2": 1.1}]
    for comp in invalid_compositions:
        with pytest.raises(ValueError):
            GasMixture("Zero/Negative mole fraction", comp, NASA_THERMAL_DATA)


def test_invalid_yaml_data():
    for invalid_data in INVALID_DATA_LIST:
        with pytest.raises(ValueError):
            GasMixture("Invalid Data", {"O2": 1.0}, invalid_data)


def test_missing_yaml_data():
    with pytest.raises(ValueError):
        GasMixture("Invalid Data", {"O2": 1.0})


def test_empty_composition():
    with pytest.raises(ValueError):
        GasMixture("Empty composition", {}, NASA_THERMAL_DATA)

    with pytest.raises(ValueError):
        GasMixture("Empty composition", thermo_data_yaml=NASA_THERMAL_DATA)


def test_empty_name():
    with pytest.raises(ValueError):
        GasMixture(
            mole_fraction_composition=STANDARD_AIR_COMPOSITION,
            thermo_data_yaml=NASA_THERMAL_DATA,
        )


def test_single_species_mixture():
    mixture = GasMixture("Pure O2", {"O2": 1.0}, NASA_THERMAL_DATA)
    assert isinstance(mixture, GasMixture)
    assert mixture.name == "Pure O2"
    assert mixture.species == {"O2"}


def test_mole_fraction_tolerance():
    valid_composition = {"N2": 0.78, "O2": 0.22 + 1e-5}  # total: 1.00001
    mixture = GasMixture(
        "Tolerant Air",
        valid_composition,
        NASA_THERMAL_DATA,
        admissible_mole_frac_error=1e-4,
    )
    assert isinstance(mixture, GasMixture)

    invalid_composition = {"N2": 0.78, "O2": 0.22 + 1e-3}  # total: 1.00022
    with pytest.raises(ValueError):
        GasMixture(
            "Too far",
            invalid_composition,
            NASA_THERMAL_DATA,
            admissible_mole_frac_error=1e-4,
        )


def test_species_case_sensitivity():
    wrong_case_composition = {"o2": 0.21, "n2": 0.79}
    with pytest.raises(ValueError):
        GasMixture("Case mismatch", wrong_case_composition, NASA_THERMAL_DATA)
