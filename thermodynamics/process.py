from enum import Enum
from typing import Optional, TypedDict, Dict, Set, Literal
from thermodynamics.state import GasState


class GasProcessType(Enum):
    ISENTROPIC = 0
    ISOTHERMAL = 1
    ISOCHORIC = 2
    ISOBARIC = 3


class GasProcessConstraint(TypedDict):
    pressure: Optional[float]
    volume: Optional[float]
    mass: Optional[float]
    temperature: Optional[float]


def process_constraint(
    pressure: Optional[float] = None,
    volume: Optional[float] = None,
    mass: Optional[float] = None,
    temperature: Optional[float] = None,
) -> GasProcessConstraint:

    constraint: GasProcessConstraint = {
        "pressure": pressure,
        "volume": volume,
        "mass": mass,
        "temperature": temperature,
    }
    num_none: int = len([val for val in constraint.values() if val is None])

    if num_none != 2:
        raise ValueError(
            "Exactly 2 of the four state variables: Pressure, Volume, Mass and Temperature, should be constrained. "
            "Except for special isentropic case, where, to constrain pressure and temperature, there should be one additional constraint on either mass or volume "
        )
    return constraint


class GasProcess:
    def __init__(
        self,
        process_type: GasProcessType,
        initial_state: GasState,
        final_state_constraints: Optional[GasProcessConstraint] = None,
    ):
        self.__process_type: GasProcessType = process_type

        self.__initial_state: GasState = initial_state
        self.__final_state: Optional[GasState] = None
        if final_state_constraints is None:
            raise ValueError(
                "Constraints for the process are required to calculate the final state"
            )
        else:
            self.__final_state_constraints: GasProcessConstraint = (
                final_state_constraints
            )
            self.__calc_final_ideal_state()

    def __validate_constraints(self):
        allowed_constraints: Dict[GasProcessType, Set[str]] = {
            GasProcessType.ISOBARIC: {"temperature", "volume", "mass"},
            GasProcessType.ISOCHORIC: {"pressure", "temperature", "mass"},
            GasProcessType.ISOTHERMAL: {"pressure", "volume", "mass"},
            GasProcessType.ISENTROPIC: {"pressure", "volume", "temperature", "mass"},
        }

        provided_constraints = {
            k for k, v in self.__final_state_constraints.items() if v is not None
        }
        invalid = provided_constraints - allowed_constraints[self.__process_type]

        if invalid:
            raise ValueError(
                f"{self.__process_type.name.capitalize()} process does not allow constraints on: "
                f"{', '.join(invalid)}"
            )

    def __calc_final_ideal_state(self):
        self.__validate_constraints()
        P1 = self.__initial_state.pressure
        V1 = self.__initial_state.volume
        M1 = self.__initial_state.mass
        T1 = self.__initial_state.temperature
        constraint_dict: GasProcessConstraint = self.__final_state_constraints
        P2 = constraint_dict.get("pressure")
        V2 = constraint_dict.get("volume")
        M2 = constraint_dict.get("mass")
        T2 = constraint_dict.get("temperature")

        if self.__process_type == GasProcessType.ISOBARIC:
            P2 = P1
            V2, M2, T2 = self.__isobaric_final_state(V1, M1, T1, V2, M2, T2)
        elif self.__process_type == GasProcessType.ISOCHORIC:
            V2 = V1
            P2, M2, T2 = self.__isochoric_final_state(P1, M1, T1, P2, M2, T2)
        elif self.__process_type == GasProcessType.ISOTHERMAL:
            T2 = T1
            P2, M2, V2 = self.__isothermal_final_state(P1, M1, V1, P2, M2, V2)
        elif self.__process_type == GasProcessType.ISENTROPIC:
            P2, V2, M2, T2 = self.__isentropic_final_state(
                P1,
                M1,
                V1,
                T1,
                P2,
                M2,
                V2,
                T2,
            )
        else:
            raise ValueError(f"{self.__process_type} process type not supported")

        self.__final_state = GasState(self.__initial_state.gas_mixture, M2, V2, P2, T2)

    def __isobaric_final_state(
        self,
        V1: float,
        M1: float,
        T1: float,
        V2: Optional[float],
        M2: Optional[float],
        T2: Optional[float],
    ):
        if len([none_val for none_val in [V2, M2, T2] if none_val is None]) != 1:
            raise ValueError(
                f"Final state information: Volume = {V2}, Mass = {M2}, Temperature = {T2}, is not appropriate for the calculation of the final isobaric state. There should be only 1 unknown variable"
            )

        def __T_relation():
            if M2 == 0 or V1 == 0 or M2 is None or V1 is None:
                raise ValueError(
                    f"For valid final state temperature calculation, final mass = {M2} or initial volume = {V1} cannot be 0"
                )
            return (M1 / M2) * T1 * (V2 / V1)

        def __M_relation():
            if V1 == 0 or T2 == 0 or V1 is None or T2 is None:
                raise ValueError(
                    f"For valid final state mass calculation, final temperature = {T2} or initial volume = {V1} cannot be 0"
                )
            return (V2 / V1) * M1 * (T1 / T2)

        def __V_relation():
            if M1 == 0 or T1 == 0 or M1 is None or T1 is None:
                raise ValueError(
                    f"For valid final state volume calculation, initial mass = {M1} or initial temperature = {T1} cannot be 0"
                )
            return (M2 / M1) * V1 * (T2 / T1)

        if T2 is None:
            T2: float = __T_relation()
        elif M2 is None:
            M2: float = __M_relation()
        elif V2 is None:
            V2: float = __V_relation()

        return V2, M2, T2

    def __isochoric_final_state(
        self,
        P1: float,
        M1: float,
        T1: float,
        P2: Optional[float],
        M2: Optional[float],
        T2: Optional[float],
    ):
        if len([none_val for none_val in [P2, M2, T2] if none_val is None]) != 1:
            raise ValueError(
                f"Final state information: Pressure = {P2}, Mass = {M2}, Temperature = {T2}, is not appropriate for the calculation of the final isochoric state. There should be only 1 unknown variable"
            )

        def __T_relation():
            if M2 == 0 or P1 == 0 or M2 is None or P1 is None:
                raise ValueError(
                    f"For valid final state temperature calculation, final mass = {M2} or initial pressure = {P1} cannot be 0"
                )
            return (M1 / M2) * T1 * (P2 / P1)

        def __M_relation():
            if T2 == 0 or P1 == 0 or T2 is None or P1 is None:
                raise ValueError(
                    f"For valid final state mass calculation, final temperature = {T2} or initial pressure = {P1} cannot be 0"
                )
            return (P2 / P1) * M1 * (T1 / T2)

        def __P_relation():
            if M1 == 0 or T1 == 0 or M1 is None or T1 is None:
                raise ValueError(
                    f"For valid final state pressure calculation, initial mass = {M1} or initial temperature = {T1} cannot be 0"
                )
            return (M2 / M1) * P1 * (T2 / T1)

        if T2 is None:
            T2: float = __T_relation()
        elif M2 is None:
            M2: float = __M_relation()
        elif P2 is None:
            P2: float = __P_relation()

        return P2, M2, T2

    def __isothermal_final_state(
        self,
        P1: float,
        M1: float,
        V1: float,
        P2: Optional[float],
        M2: Optional[float],
        V2: Optional[float],
    ):
        if len([none_val for none_val in [P2, M2, V2] if none_val is None]) != 1:
            raise ValueError(
                f"Final state information: Pressure = {P2}, Mass = {M2}, Volume = {V2}, is not appropriate for the calculation of the final isothermal state. There should be only 1 unknown variable"
            )

        def __M_relation():
            if P1 == 0 or V1 == 0 or P1 is None or V1 is None:
                raise ValueError(
                    f"For valid final state mass calculation, initial pressure = {P1} or initial volume = {V1} cannot be 0"
                )
            return (V2 / V1) * M1 * (P2 / P1)

        def __V_relation():
            if P2 == 0 or M1 == 0 or P2 is None or M1 is None:
                raise ValueError(
                    f"For valid final state volume calculation, final pressure = {P2} or initial mass = {M1} cannot be 0"
                )
            return (P1 / P2) * V1 * (M2 / M1)

        def __P_relation():
            if V2 == 0 or M1 == 0 or V2 is None or M1 is None:
                raise ValueError(
                    f"For valid final state pressure calculation, final volume = {V2} or initial mass = {M1} cannot be 0"
                )
            return (V1 / V2) * P1 * (M2 / M1)

        if M2 is None:
            M2: float = __M_relation()
        elif V2 is None:
            V2: float = __V_relation()
        elif P2 is None:
            P2: float = __P_relation()
        return P2, M2, V2

    def __isentropic_final_state(
        self,
        P1: float,
        M1: float,
        V1: float,
        T1: float,
        P2: Optional[float],
        M2: Optional[float],
        V2: Optional[float],
        T2: Optional[float],
    ):
        y = self.__initial_state.gamma
        if y <= 0:
            raise ValueError(
                f"For an isentropic process, gamma = {y} must be a positive non-zero value"
            )
        if len([none_val for none_val in [P2, M2, V2, T2] if none_val is None]) not in [
            2,
            3,
        ]:
            raise ValueError(
                f"Final state information: Pressure = {P2}, Mass = {M2}, Volume = {V2}, Temperature = {T2},"
                " is not appropriate for the calculation of the final isentropic state. "
                "There should be only 2 constraints except for a special isentropic case, where, to constrain pressure and temperature, "
                "there should be one additional constraint on either mass or volume"
            )

        def __V_M_P_relation():
            if P2 == 0 or P2 is None or M1 == 0 or M1 is None:
                raise ValueError(
                    f"For valid final state volume calculation, final pressure = {P2} or initial mass = {M1} cannot be 0"
                )
            return ((P1 / P2) ** (1 / y)) * (M2 / M1) * V1

        def __T_V_M_relation():
            if V2 == 0 or V2 is None or M1 == 0 or M1 is None:
                raise ValueError(
                    f"For valid final state temperature calculation, final volume = {V2} or initial mass = {M1} cannot be 0"
                )
            return T1 * ((V1 / V2) ** (y - 1)) * ((M2 / M1) ** (y - 1))

        def __M_P_V_relation():
            if P1 == 0 or P1 is None or V1 == 0 or V1 is None:
                raise ValueError(
                    f"For valid final state mass calculation, initial pressure = {P1} or initial volume = {V1} cannot be 0"
                )
            return ((P2 / P1) ** (1 / y)) * (V2 / V1) * M1

        def __P_M_V_relation():
            if M1 == 0 or M1 is None or V2 == 0 or V2 is None:
                raise ValueError(
                    f"For valid final state pressure calculation, initial mass = {M1} or final volume = {V2} cannot be 0"
                )
            return ((M2 / M1) ** y) * ((V1 / V2) ** y) * P1

        def __M_V_T_relation():
            if V1 == 0 or V1 is None or T1 == 0 or T1 is None:
                raise ValueError(
                    f"For valid final state mass calculation, initial volume = {V1} or initial temperature = {T1} cannot be 0"
                )
            return (V2 / V1) * ((T2 / T1) ** (1 / (y - 1))) * M1

        def __V_T_M_relation():
            if T2 == 0 or T2 is None or M1 == 0 or M1 is None:
                raise ValueError(
                    f"For valid final state volume calculation, final temperature = {T2} or initial mass = {M1} cannot be 0"
                )
            return (T1 / T2) ** (1 / (y - 1)) * (M2 / M1) * V1

        if V2 is None and T2 is None and P2 is not None and M2 is not None:
            V2 = __V_M_P_relation()
            T2 = __T_V_M_relation()
        elif M2 is None and T2 in None and P2 is not None and V2 is not None:
            M2 = __M_P_V_relation()
            T2 = __T_V_M_relation()
        elif P2 is None and T2 is None and M2 is not None and V2 is not None:
            P2 = __P_M_V_relation()
            T2 = __T_V_M_relation()
        elif M2 is None and P2 is None and V2 is not None and T2 is not None:
            M2 = __M_V_T_relation()
            P2 = __P_M_V_relation()
        elif V2 is None and P2 is None and M2 is not None and T2 is not None:
            V2 = __V_T_M_relation()
            P2 = __P_M_V_relation()
        elif V2 is None and M2 is not None and P2 is not None and T2 is not None:
            V2 = __V_M_P_relation()
        elif M2 is None and V2 is not None and P2 is not None and T2 is not None:
            M2 = __M_V_T_relation()
        else:
            raise ValueError(
                f"Invalid Constraints: Pressure = {P2}, Mass = {M2}, Volume = {V2}, Temperature = {T2},"
                " for the calculation of the final isentropic state. "
                "There should be only 2 constraints except for a special isentropic case, where, to constrain pressure and temperature, "
                "there should be one additional constraint on either mass or volume"
            )

        return P2, V2, M2, T2

    def __traverse_to_final_state_variable(
        self,
        state_variable: Literal["pressure", "volume", "mass", "temperature"],
        process_type: GasProcessType,
        num_steps: int,
        gamma_update: Literal["constant", "last", "mean", "mean-global"],
    ):
        raise NotImplementedError("Process interpolation is not implemented yet.")

    @property
    def initial_state(self):
        return self.__initial_state

    @property
    def final_state(self):
        return self.__final_state
