import math
from typing import Dict, Set, List, Tuple
from thermodynamics.phase import GasMixture
from cache.CONSTANTS import UNIVERSAL_GAS_CONSTANT_SI as R_SI


class GasState:
    def __init__(
        self,
        gas_mixture: GasMixture,
        mass_kg: float = None,
        volume_SI: float = None,
        pressure_Pa: float = None,
        temperature_K: float = None,
    ):
        """GasState is an abstraction of a thermodynamic state of a gas mixture.
        Only 3 of the 4 state variables: Mass, Volume, Pressure, Temperature need to be provided to
        completely define the state.
        However, if all 4 are provided, they are validated using the equation of state, failing which a ValueError is raised

        Args:
            gas_mixture (GasMixture): Contains the thermal data of the mixture species
            mass_kg (float, optional): Mass of the gas mixture. Defaults to None.
            volume_SI (float, optional): Volume of the gas mixture. Defaults to None.
            pressure_Pa (float, optional): Pressure of the gas mixture. Defaults to None.
            temperature_K (float, optional): Temperature of the gas. Defaults to None.
        """
        # Gas mixture to be evaluated
        self.__gas_mixture: GasMixture = gas_mixture

        # Molar mass of the mixture in kg
        self.__mol_mass_kg: float = None

        # Mass of the mixture in kg
        self.__mass_kg: float = mass_kg

        self.__volume_SI: float = volume_SI

        # Gas Pressure
        self.__pressure_Pa: float = pressure_Pa

        # Gas temperature
        self.__temperature_K: float = temperature_K

        self.__calc_molar_mass()

        # Specfic gas constant of the gas species
        self.__R_spec = R_SI / self.__mol_mass_kg

        self.__eqn_of_state()

        self.__validate_eqn_of_state()

        # Isobaric molar specific heat
        self.__cp_SI: float = None

        # Isochoric molar specific heat
        self.__cv_SI: float = None

        # Molar enthalpy
        self.__enthalpy_SI: float = None

        # Molar entropy
        self.__entropy_SI: float = None

        # Ratio of specific heats
        self.__gamma: float = None

        # Calculates all the thermodynamic properties for a given state
        self.__calc_thermodynamic_properties()

    def __calc_molar_mass(self):
        """Calculates molar mass of the entire mixture"""
        molar_mass_total = 0.0
        for species, mole_frac in self.__gas_mixture.mole_fraction_composition.items():
            species_data = self.__gas_mixture.species_data[species]
            species_mol_mass_kg = species_data["molar-mass"]
            molar_mass = mole_frac * species_mol_mass_kg
            molar_mass_total += molar_mass

        self.__mol_mass_kg = molar_mass_total

    def __eqn_of_state(self):
        """Equation of state `Pressure * Volume = mass * specific gas constant * temperature`

        Raises:
            ValueError: If out of 4 state variables mass, pressure, volume, temperature, 3 are not provided.
        """
        # Pressure is unknown
        if (
            self.__pressure_Pa is None
            and self.__volume_SI is not None
            and self.__mass_kg is not None
            and self.__temperature_K is not None
        ):
            self.__pressure_Pa = (
                self.__mass_kg * self.__R_spec * self.__temperature_K
            ) / self.__volume_SI
        # Volume is unknown
        elif (
            self.__pressure_Pa is not None
            and self.__volume_SI is None
            and self.__mass_kg is not None
            and self.__temperature_K is not None
        ):
            self.__volume_SI = (
                self.__mass_kg * self.__R_spec * self.__temperature_K
            ) / self.__pressure_Pa

        # Mass is unknown
        elif (
            self.__pressure_Pa is not None
            and self.__volume_SI is not None
            and self.__mass_kg is None
            and self.__temperature_K is not None
        ):
            self.__mass_kg = (self.__pressure_Pa * self.__volume_SI) / (
                self.__R_spec * self.__temperature_K
            )
        # Temperature is unknown
        elif (
            self.__pressure_Pa is not None
            and self.__volume_SI is not None
            and self.__mass_kg is not None
            and self.__temperature_K is None
        ):
            self.__temperature_K = (self.__pressure_Pa * self.__volume_SI) / (
                self.__R_spec * self.__mass_kg
            )
        elif (
            self.__pressure_Pa is not None
            and self.__volume_SI is not None
            and self.__mass_kg is not None
            and self.__temperature_K is not None
        ):
            self.__validate_eqn_of_state()
        else:
            raise ValueError(
                "From mass, volume, pressure and temperature, atleast 3 must be be provided for a valid calculation of the state"
            )

    def __validate_eqn_of_state(self, admissible_error: float = 1e-6):
        """Validate if all the state variables are actually valid

        Args:
            admissible_error (float, optional): Maximum admissible error between the lhs and rhs of equation of state.. Defaults to 1e-6.

        Raises:
            ValueError: If state variables, mass, volume, pressure and temperature are not valid
        """
        difference = (self.__pressure_Pa * self.__volume_SI) - (
            self.__mass_kg * self.__R_spec * self.__temperature_K
        )
        validation_criteria: bool = abs(difference) <= admissible_error

        if not validation_criteria:
            raise ValueError(
                f"State with pressure {self.__pressure_Pa}Pa, volume {self.__volume_SI},mass {self.__mass_kg} and temperature {self.__temperature_K} is not valid thermodynamically. L.H.S Eqn of State - R.H.S of Eqn of state: {abs(difference)}. Validation criteria: {admissible_error}"
            )

    def __thermo_nasa_polynomials(
        self, coefficients: List[float]
    ) -> Tuple[float, float, float]:
        """NASA polynomials for the thermodynamic quantities of the species.

        Args:
            coefficients (List[float]): 7 NASA coefficicents of the species in a specified temperature range

        Returns:
            (float, float, float): cp / R, h / RT, s / R
        """
        t = self.__temperature_K
        a1, a2, a3, a4, a5, a6, a7 = coefficients
        cp_R = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4
        h_RT = (
            a1
            + (a2 * t / 2)
            + (a3 * t**2 / 3)
            + (a4 * t**3 / 4)
            + (a5 * t**4 / 5)
            + (a6 / t)
        )

        s_R = (
            a1 * math.log(t)
            + a2 * t
            + (a3 / 2) * t**2
            + (a4 / 3) * t**3
            + (a5 / 4) * t**4
            + a7
        )
        return cp_R, h_RT, s_R

    def __calc_thermodynamic_properties(self):
        """Calculates all the relevant thermodynamic properties of the gas mixture

        Raises:
            ValueError: If the temperature of the state lies beyond the scope of the temperature ranges in the parsed yaml file
        """
        cp_total = 0.0
        enthalpy_total = 0.0
        entropy_total = 0.0
        gas_mixture_species: Set[str] = self.__gas_mixture.species

        for species in gas_mixture_species:
            species_data: Dict[str, Dict] = self.__gas_mixture.species_data[species]
            temp_range = species_data["temperature-ranges"]
            coeff_array = species_data["data"]
            temp_range_len = len(temp_range)

            if temp_range_len == 3:
                t_lo, t_mid, t_hi = temp_range

                if t_lo <= self.__temperature_K <= t_mid:
                    coeffs = coeff_array[0]

                elif t_mid < self.__temperature_K <= t_hi:
                    coeffs = coeff_array[1]

                else:
                    raise ValueError(
                        f"[{species}] T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            else:
                t_lo, t_hi = temp_range

                if t_lo <= self.__temperature_K <= t_hi:
                    coeffs = coeff_array[0]
                else:
                    raise ValueError(
                        f"[{species}] T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            cp_R, h_RT, s_R = self.__thermo_nasa_polynomials(coeffs)
            mole_frac = self.__gas_mixture.mole_fraction_composition[species]
            cp = cp_R * mole_frac
            enthalpy = h_RT * mole_frac
            entropy = s_R * mole_frac

            cp_total += cp
            enthalpy_total += enthalpy
            entropy_total += entropy

        self.__cp_SI = cp_total * R_SI
        self.__enthalpy_SI = enthalpy_total * R_SI * self.__temperature_K
        self.__entropy_SI = entropy_total * R_SI
        self.__cv_SI = self.__cp_SI - R_SI
        self.__gamma = self.__cp_SI / self.__cv_SI

    def update_P_V_T(self, pressure_Pa: float, volume_SI: float, temperature_K: float):
        """Updates Pressure, Volume and Temperature of the state

        Args:
            pressure_Pa (float): Pressure in Pascal
            volume_SI (float): Volume in cubic meters
            temperature_K (float): Temperature in Kelvin
        """
        self.__mass_kg = None
        self.__pressure_Pa = pressure_Pa
        self.__volume_SI = volume_SI
        self.__temperature_K = temperature_K
        self.__eqn_of_state()
        self.__validate_eqn_of_state()
        self.__calc_thermodynamic_properties()

    def update_P_m_T(self, pressure_Pa: float, mass_kg: float, temperature_K: float):
        """Updates Pressure, Mass and Temperature of the state

        Args:
            pressure_Pa (float): Pressure in Pascal
            mass_kg (float): Mass of the gas mixture in kilograms
            temperature_K (float): Temperature in Kelvin
        """
        self.__volume_SI = None
        self.__pressure_Pa = pressure_Pa
        self.__mass_kg = mass_kg
        self.__temperature_K = temperature_K
        self.__eqn_of_state()
        self.__validate_eqn_of_state()
        self.__calc_thermodynamic_properties()

    def update_V_m_T(self, volume_SI: float, mass_kg: float, temperature_K: float):
        """Updates Volume, Mass and Temperature of the state

        Args:
            volume_SI (float): Volume in cubic meters
            mass_kg (float): Mass of the gas mixture in kilograms
            temperature_K (float): Temperature in Kelvin
        """
        self.__pressure_Pa = None
        self.__mass_kg = mass_kg
        self.__volume_SI = volume_SI
        self.__temperature_K = temperature_K
        self.__eqn_of_state()
        self.__validate_eqn_of_state()
        self.__calc_thermodynamic_properties()

    def update_P_V_m(self, pressure_Pa: float, volume_SI: float, mass_kg: float):
        """Updates Pressure, Volume and Mass of the state

        Args:
            pressure_Pa (float): Pressure in Pascal
            volume_SI (float): Volume in cubic meters
            mass_kg (float): Mass of the gas mixture in kilograms
        """
        self.__temperature_K = None
        self.__pressure_Pa = pressure_Pa
        self.__volume_SI = volume_SI
        self.__mass_kg = mass_kg
        self.__eqn_of_state()
        self.__validate_eqn_of_state()
        self.__calc_thermodynamic_properties()

    @property
    def cp(self) -> float:
        return self.__cp_SI / self.molar_mass

    @property
    def cv(self) -> float:
        return self.__cv_SI / self.molar_mass

    @property
    def enthalpy(self) -> float:
        return self.__enthalpy_SI / self.molar_mass

    @property
    def entropy(self) -> float:
        return self.__entropy_SI / self.molar_mass

    @property
    def gamma(self) -> float:
        return self.__gamma

    @property
    def molar_mass(self) -> float:
        return self.__mol_mass_kg

    @property
    def pressure(self) -> float:
        return self.__pressure_Pa

    @property
    def density(self) -> float:
        return self.__mass_kg / self.__volume_SI

    @property
    def mass(self) -> float:
        return self.__mass_kg

    @property
    def volume(self) -> float:
        return self.__volume_SI

    @property
    def R(self) -> float:
        return self.__R_spec

    @property
    def R_uni(self) -> float:
        return R_SI

    @property
    def temperature(self) -> float:
        return self.__temperature_K

    @property
    def gas_mixture(self) -> GasMixture:
        return self.__gas_mixture

    def __str__(self):
        lines = [self.__gas_mixture.__str__()]
        lines.append(f"\nTHERMODYNAMIC DATA FOR '{self.__gas_mixture.name}':\n")
        lines.append(f"Molar Mass (kg): {self.molar_mass}")
        lines.append(f"Mass (kg): {self.mass}")
        lines.append(f"Pressure (Pa): {self.pressure}")
        lines.append(f"Volume (cub. m): {self.volume}")
        lines.append(f"Density (kg / cub. m): {self.density}")
        lines.append(f"Temperature (K): {self.temperature}")
        lines.append(f"Universal Gas Constant (J / mol K): {self.R_uni}")
        lines.append(f"Specific Gas Constant (J / kg K): {self.R}")
        lines.append(f"Isobaric specific heat  (J / kg K): {self.cp}")
        lines.append(f"Isochoric specific heat (J / kg K): {self.cv}")
        lines.append(f"Specific Enthalpy (J / kg): {self.enthalpy}")
        lines.append(f"Specific Entropy (J / kg K): {self.entropy}")
        lines.append(f"Ratio of specific heats: {self.gamma}")

        return "\n".join(lines)
