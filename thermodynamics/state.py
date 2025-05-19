import math
from typing import Dict, Set, List, Tuple, Optional, Any
from thermodynamics.phase import GasMixture
from cache.CONSTANTS import UNIVERSAL_GAS_CONSTANT_SI as R_SI


class GasState:
    def __init__(
        self,
        name: Optional[str] = None,
        gas_mixture: GasMixture = None,
        pressure_Pa: Optional[float] = None,
        density_SI: Optional[float] = None,
        temperature_K: Optional[float] = None,
    ):
        """GasState is an abstraction of a thermodynamic state of a gas mixture.
        Only 2 of the 3 state variables: Pressure, Density and Temperature need to be provided to
        completely define the state.
        However, if all 3 are provided, they are validated using the equation of state, failing which an error is raised

        Args:
            name (str): Name of the thermodynamic state.
            gas_mixture (GasMixture): Contains the thermal data of the mixture species
            pressure_Pa (float, optional): Pressure of the gas mixture. Defaults to None.
            density_SI (float, optional): Density of the gas mixture. Defaults to None.
            temperature_K (float, optional): Temperature of the gas. Defaults to None.

        Raises:
            ValueError: If name of the state is not provided.
            ValueError: If out of 3 state variables pressure, density and temperature, 2 are not provided.
            ValueError: If state variables pressure, density and temperature are not valid.
            ValueError: If the temperature of the state lies beyond the scope of the temperature ranges in the parsed yaml file.
            ValueError: If any of the state variables is negative or zero.
        """
        if name is None:
            raise ValueError("Name must be provided for the state")
        if gas_mixture == None:
            raise ValueError("Gas mixture must be provided for the state")

        self.__name = name

        # Gas mixture to be evaluated
        self.__gas_mixture: GasMixture = gas_mixture

        # Molar mass of the mixture in kg
        self.__mol_mass_kg: Optional[float] = None

        # Mass of the mixture in kg
        self.__density_SI = density_SI

        # Gas Pressure
        self.__pressure_Pa: Optional[float] = pressure_Pa

        # Gas temperature
        self.__temperature_K: Optional[float] = temperature_K

        self.__calc_molar_mass()

        # Specfic gas constant of the gas species
        self.__R_spec: float = R_SI / self.__mol_mass_kg

        self.__eqn_of_state()

        # Isobaric molar specific heat
        self.__mol_cp_SI: Optional[float] = None

        # Isochoric molar specific heat
        self.__mol_cv_SI: Optional[float] = None

        # Molar enthalpy
        self.__mol_enthalpy_SI: Optional[float] = None

        # Molar entropy
        self.__mol_entropy_SI: Optional[float] = None

        # Ratio of specific heats
        self.__gamma: Optional[float] = None

        # Calculates all the thermodynamic properties for a given state
        self.__calc_thermodynamic_properties()

    def __calc_molar_mass(self):
        """Calculates molar mass of the entire mixture"""
        molar_mass_total: float = 0.0
        for species, mole_frac in self.__gas_mixture.mole_fraction_composition.items():
            species_data = self.__gas_mixture.species_data[species]
            species_mol_mass_kg = species_data["molar-mass"]
            molar_mass = mole_frac * species_mol_mass_kg
            molar_mass_total += molar_mass

        self.__mol_mass_kg = molar_mass_total

    def __eqn_of_state(self):
        """Equation of state `pressure = density * specific gas constant * temperature`

        Raises:
            ValueError: If out of 3 state variables pressure, density and temperature, 2 are not provided.
            ValueError: If any of the state variables is negative or zero.
        """
        if any(
            x is not None and x <= 0
            for x in [self.__pressure_Pa, self.__density_SI, self.__temperature_K]
        ):
            raise ValueError("None of the state variables can be negative or zero")
        # Pressure is unknown
        if (
            self.__pressure_Pa is None
            and self.__density_SI is not None
            and self.__temperature_K is not None
        ):
            self.__pressure_Pa = (
                self.__density_SI * self.__R_spec * self.__temperature_K
            )
        # Density is unknown
        elif (
            self.__pressure_Pa is not None
            and self.__density_SI is None
            and self.__temperature_K is not None
        ):
            self.__density_SI = self.__pressure_Pa / (
                self.__R_spec * self.__temperature_K
            )

        # Temperature is unknown
        elif (
            self.__pressure_Pa is not None
            and self.__density_SI is not None
            and self.__temperature_K is None
        ):
            self.__temperature_K = (self.__pressure_Pa) / (
                self.__R_spec * self.__density_SI
            )
        # Entire state is provided
        elif (
            self.__pressure_Pa is not None
            and self.__density_SI is not None
            and self.__temperature_K is not None
        ):
            self.__validate_eqn_of_state()
        else:
            raise ValueError(
                "From pressure, density and temperature, atleast 2 must be be provided for a valid calculation of the state"
            )

    def __validate_eqn_of_state(self, admissible_error: float = 1e-5):
        """Validate if all the state variables are actually valid

        Args:
            admissible_error (float, optional): Maximum admissible relative error between the lhs and rhs of equation of state.. Defaults to 1e-6.

        Raises:
            ValueError: If state variables pressure, density and temperature are not valid
        """
        difference: Optional[float] = (self.__pressure_Pa) - (
            self.__density_SI * self.__R_spec * self.__temperature_K
        )
        validation_criteria: bool = (
            abs(difference / self.__pressure_Pa) <= admissible_error
        )

        if not validation_criteria:
            raise ValueError(
                f"State with pressure {self.__pressure_Pa}Pa, density {self.__density_SI} and temperature {self.__temperature_K} is not valid thermodynamically. Relative Error: {abs((difference)/ self.__pressure_Pa)}. Validation criteria: {admissible_error}"
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
        t: float = self.__temperature_K
        a1, a2, a3, a4, a5, a6, a7 = coefficients
        cp_R: float = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4
        h_RT: float = (
            a1
            + (a2 * t / 2)
            + (a3 * t**2 / 3)
            + (a4 * t**3 / 4)
            + (a5 * t**4 / 5)
            + (a6 / t)
        )

        s_R: float = (
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
        cp_total: float = 0.0
        enthalpy_total: float = 0.0
        entropy_total: float = 0.0
        gas_mixture_species: Set[str] = self.__gas_mixture.species

        for species in gas_mixture_species:
            species_data: Dict[str, Any] = self.__gas_mixture.species_data[species]
            temp_range: List[float] = species_data["temperature-ranges"]
            coeff_array: List[List[float]] = species_data["data"]
            temp_range_len: int = len(temp_range)

            if temp_range_len == 3:
                t_lo, t_mid, t_hi = temp_range

                if t_lo <= self.__temperature_K <= t_mid:
                    coeffs = coeff_array[0]

                elif t_mid <= self.__temperature_K <= t_hi:
                    coeffs = coeff_array[1]

                else:
                    raise ValueError(
                        f"For species {species}: T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            else:
                t_lo, t_hi = temp_range

                if t_lo <= self.__temperature_K <= t_hi:
                    coeffs = coeff_array[0]
                else:
                    raise ValueError(
                        f"For species {species}: T={self.__temperature_K}K not in range ({t_lo}-{t_hi})K"
                    )

            cp_R, h_RT, s_R = self.__thermo_nasa_polynomials(coeffs)
            mole_frac = self.__gas_mixture.mole_fraction_composition[species]
            cp = cp_R * mole_frac
            enthalpy = h_RT * mole_frac
            entropy = s_R * mole_frac

            cp_total += cp
            enthalpy_total += enthalpy
            entropy_total += entropy

        self.__mol_cp_SI = cp_total * R_SI
        self.__mol_enthalpy_SI = enthalpy_total * R_SI * self.__temperature_K
        self.__mol_entropy_SI = entropy_total * R_SI
        self.__mol_cv_SI = self.__mol_cp_SI - R_SI
        self.__gamma = self.__mol_cp_SI / self.__mol_cv_SI

    def update_state(
        self,
        name: Optional[str] = None,
        pressure_Pa: Optional[float] = None,
        density_SI: Optional[float] = None,
        temperature_K: Optional[float] = None,
    ):
        """Updates the state for the two provided state variables.

        Args:
            name (Optional[str], optional): Updated name. Defaults to existing state name.
            pressure_Pa (Optional[float], optional): Updated pressure. Defaults to None.
            density_SI (Optional[float], optional): Updated density. Defaults to None.
            temperature_K (Optional[float], optional): Updated temperature. Defaults to None.

        Raises:
            ValueError: If no updates are provided
            ValueError: If only a single state variable is provided as an update
        """
        if name == None:
            name = self.__name
        if (
            pressure_Pa is None and density_SI is None and temperature_K is None
        ):  # All updates are missing
            raise ValueError("No inputs provided for state update")
        elif (
            (pressure_Pa is None and density_SI is None)
            or (pressure_Pa is None and temperature_K is None)
            or (density_SI is None and temperature_K is None)
        ):  # Only one update is provided
            raise ValueError(
                "2 state variables out of the 3: pressure, density and temperature, need to be updated"
            )
        self.__name = name
        self.__pressure_Pa = pressure_Pa
        self.__density_SI = density_SI
        self.__temperature_K = temperature_K
        self.__eqn_of_state()
        self.__calc_thermodynamic_properties()

    @property
    def cp_spec(self) -> float:
        return float(self.__mol_cp_SI / self.molar_mass)

    @property
    def cp_mol(self) -> float:
        return float(self.__mol_cp_SI)

    @property
    def cv_spec(self) -> float:
        return float(self.__mol_cv_SI / self.molar_mass)

    @property
    def cv_mol(self) -> float:
        return float(self.__mol_cv_SI)

    @property
    def enthalpy_spec(self) -> float:
        return float(self.__mol_enthalpy_SI / self.molar_mass)

    @property
    def enthalpy_mol(self) -> float:
        return float(self.__mol_enthalpy_SI)

    @property
    def state_name(self) -> str:
        return self.__name

    @property
    def entropy_spec(self) -> float:
        return float(self.__mol_entropy_SI / self.molar_mass)

    @property
    def entropy_mol(self) -> float:
        return float(self.__mol_entropy_SI)

    @property
    def gamma(self) -> float:
        return float(self.__gamma)

    @property
    def molar_mass(self) -> float:
        return float(self.__mol_mass_kg)

    @property
    def pressure(self) -> float:
        return float(self.__pressure_Pa)

    @property
    def density(self) -> float:
        return float(self.__density_SI)

    @property
    def R_spec(self) -> float:
        return float(self.__R_spec)

    @property
    def R(self) -> float:
        return float(R_SI)

    @property
    def temperature(self) -> float:
        return float(self.__temperature_K)

    @property
    def gas_mixture(self) -> GasMixture:
        return self.__gas_mixture

    def __str__(self):
        # lines = [self.__gas_mixture.__str__()]
        lines = []
        lines.append(
            f"\nTHERMODYNAMIC DATA OF '{self.__gas_mixture.mixture_name}' AT STATE: {self.state_name}\n"
        )
        lines.append(f"Molar Mass (kg): {self.molar_mass}")
        lines.append(f"Pressure (Pa): {self.pressure}")
        lines.append(f"Density (kg / cub. m): {self.density}")
        lines.append(f"Temperature (K): {self.temperature}")
        lines.append(f"Universal Gas Constant (J / mol K): {self.R}")
        lines.append(f"Specific Gas Constant (J / kg K): {self.R_spec}")
        lines.append(f"Specific Cp (J / kg K): {self.cp_spec}")
        lines.append(f"Molar Cp (J / mol K): {self.cp_mol}")
        lines.append(f"Specific Cv (J / kg K): {self.cv_spec}")
        lines.append(f"Molar Cv (J / mol K): {self.cv_mol}")
        lines.append(f"Specific Enthalpy (J / kg): {self.enthalpy_spec}")
        lines.append(f"Molar Enthalpy (J / mol): {self.enthalpy_mol}")
        lines.append(f"Specific Entropy (J / kg K): {self.entropy_spec}")
        lines.append(f"Molar Entropy (J / mol K): {self.entropy_mol}")
        lines.append(f"Ratio of specific heats: {self.gamma}")

        return "\n".join(lines)
