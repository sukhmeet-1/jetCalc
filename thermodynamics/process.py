from enum import Enum
from typing import Optional, TypedDict
from thermodynamics.state import GasState


class GasProcessType(Enum):
    ADIABATIC = 0
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
            "Exactly 2 of the four state variables: Pressure, Volume, Mass and Temperature, should be constrained"
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
            self.__calc_final_state()

    def __calc_final_state(self):
        allowed_constraints = {
            GasProcessType.ISOBARIC: {"temperature", "volume", "mass"},
            GasProcessType.ISOCHORIC: {"pressure", "temperature", "mass"},
            GasProcessType.ISOTHERMAL: {"pressure", "volume", "mass"},
            GasProcessType.ADIABATIC: {"pressure", "volume", "temperature", "mass"},
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
        elif self.__process_type == GasProcessType.ADIABATIC:
            P2, V2, M2, T2 = self.__adiabatic_final_state(
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
        if V2 is not None and M2 is not None:
            if M2 == 0 or V1 == 0:
                raise ValueError(
                    "For valid final state temperature calculation, final mass or initial volume cannot be 0"
                )
            T2: float = (M1 / M2) * T1 * (V2 / V1)
        elif V2 is not None and T2 is not None:
            if V1 == 0 or T2 == 0:
                raise ValueError(
                    "For valid final state mass calculation, final temperature or initial volume cannot be 0"
                )
            M2: float = (V2 / V1) * M1 * (T1 / T2)
        elif T2 is not None and M2 is not None:
            if M1 == 0 or T1 == 0:
                raise ValueError(
                    "For valid final state volume calculation, initial mass or initial temperature cannot be 0"
                )
            V2: float = (M2 / M1) * V1 * (T2 / T1)
        else:
            raise ValueError(
                "Isobaric process can only constrain any 2 of Temperature, Volume and Mass."
            )

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
        if P2 is not None and M2 is not None:
            if M2 == 0 or P1 == 0:
                raise ValueError(
                    "For valid final state temperature calculation, final mass or initial pressure cannot be 0"
                )
            T2: float = (M1 / M2) * T1 * (P2 / P1)
        elif P2 is not None and T2 is not None:
            if T2 == 0 or P1 == 0:
                raise ValueError(
                    "For valid final state mass calculation, final temperature or initial pressure cannot be 0"
                )
            M2: float = (P2 / P1) * M1 * (T1 / T2)
        elif T2 is not None and M2 is not None:
            if M1 == 0 or T1 == 0:
                raise ValueError(
                    "For valid final state pressure calculation, initial mass or initial temperature cannot be 0"
                )
            P2: float = (M2 / M1) * P1 * (T2 / T1)
        else:
            raise ValueError(
                "Isochoric process can only constrain any 2 of Pressure, Temperature and Mass."
            )

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
        if P2 is not None and V2 is not None:
            if P1 == 0 or V1 == 0:
                raise ValueError(
                    "For valid final state mass calculation, initial pressure or initial volume cannot be 0"
                )
            M2: float = (V2 / V1) * M1 * (P2 / P1)
        elif P2 is not None and M2 is not None:
            if P2 == 0 or M1 == 0:
                raise ValueError(
                    "For valid final state volume calculation, final pressure or initial mass cannot be 0"
                )
            V2: float = (P1 / P2) * V1 * (M2 / M1)
        elif V2 is not None and M2 is not None:
            if V2 == 0 or M1 == 0:
                raise ValueError(
                    "For valid final state pressure calculation, final volume or initial mass cannot be 0"
                )
            P2: float = (V1 / V2) * P1 * (M2 / M1)
        else:
            raise ValueError(
                "Isothermal process can only constrain any 2 of Pressure, Volume and Mass."
            )

        return P2, M2, V2

    def __adiabatic_final_state(
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
        if y == 0 or y < 0:
            raise ValueError(
                "For an adiabatic process, gamma must have a positive non-zero value"
            )
        if P2 is not None and M2 is not None:
            if P2 == 0:
                raise ValueError(
                    "For valid final state volume calculation, final pressure cannot be 0"
                )
            V2 = ((P1 / P2) ** (1 / y)) * V1
            if M2 == 0 or V2 == 0:
                raise ValueError(
                    "For valid final state temperature calculation, final mass or final volume cannot be 0"
                )
            T2 = (M1 / M2) * T1 * (V1 / V2) ** (y - 1)
        elif P2 is not None and T2 is not None:
            if P2 == 0:
                raise ValueError(
                    "For valid final state volume calculation, final pressure cannot be 0"
                )
            V2 = ((P1 / P2) ** (1 / y)) * V1
            if T2 == 0 or V2 == 0:
                raise ValueError(
                    "For valid final state mass calculation, final temperature or final volume cannot be 0"
                )
            M2 = (T1 / T2) * M1 * ((V1 / V2) ** (y - 1))
        elif M2 is not None and T2 is not None:
            if M2 == 0 or T2 == 0:
                raise ValueError(
                    "For valid final state volume calculation, final mass or initial temperature cannot be 0"
                )
            V2 = (M1 / M2) * (T1 / T2) * (V1 ** (y - 1))
            if V2 == 0:
                raise ValueError(
                    "For valid final state pressure calculation, final volume cannot be 0"
                )
            P2 = P1 * ((V1 / V2) ** y)
        elif M2 is not None and V2 is not None:
            if V2 == 0:
                raise ValueError(
                    "For valid final state pressure calculation, final volume cannot be 0"
                )
            P2 = ((V1 / V2) ** y) * P1
            if M2 == 0 or V2 == 0:
                raise ValueError(
                    "For valid final state temperature calculation, final mass or final volume cannot be 0"
                )
            T2 = (M1 / M2) * T1 * ((V1 / V2) ** (y - 1))
        else:
            raise ValueError(
                "Adiabatic process can constrain any 2 of Pressure, Volume, Temperature and Mass."
            )

        return P2, V2, M2, T2

    def __interpolate_process(self):
        raise NotImplementedError("Process interpolation is not implemented yet.")

    @property
    def initial_state(self):
        return self.__initial_state

    @property
    def final_state(self):
        return self.__final_state
