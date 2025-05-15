import math
from typing import Literal, Dict
from thermodynamics.state import GasState

PROCESS_TYPES = ["adiabatic", "isothermal", "isochoric", "isobaric"]


class GasProcess:
    def __init__(
        self,
        process_type: Literal["adiabatic", "isothermal", "isochoric", "isobaric"],
        initial_state: GasState = None,
        process_constraints: Dict[str, float] = None,
        final_state: GasState = None,
    ):
        if process_type not in PROCESS_TYPES:
            raise ValueError(
                f"'{process_type}' process type is not supported. Supported types are {PROCESS_TYPES}"
            )
        self.__process_type: str = process_type
        if initial_state is None:
            raise ValueError("Initial state must be provided.")
        self.__initial_state: GasState = initial_state
        self.__final_state: GasState = final_state
        self.__final_state_constraints = process_constraints
        assert (
            len(self.__final_state_constraints.keys()) == 2
        ), "Process should have only 2 constraints out of the 4 state variables: Pressure, Volume, Mass, Temperature"
        if self.__final_state_constraints is not None:
            self.__calc_final_state()

    def __calc_final_state(self):
        P1 = self.__initial_state.pressure
        V1 = self.__initial_state.volume
        M1 = self.__initial_state.mass
        T1 = self.__initial_state.temperature
        constraint_dict: Dict[str, float] = self.__final_state_constraints
        P2 = constraint_dict.get("Pressure", None)
        V2 = constraint_dict.get("Volume", None)
        M2 = constraint_dict.get("Mass", None)
        T2 = constraint_dict.get("Temperature", None)

        if self.__process_type == "isobaric":
            P2 = P1
            V2, M2, T2 = self.__isobaric_final_state(V1, M1, T1, V2, M2, T2)
        elif self.__process_type == "isochoric":
            V2 = V1
            P2, M2, T2 = self.__isochoric_final_state(P1, M1, T1, P2, M2, T2)
        elif self.__process_type == "isothermal":
            T2 = T1
            P2, M2, V2 = self.__isothermal_final_state(P1, M1, V1, P2, M2, V2)
        elif self.__process_type == "adiabatic":
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
        self, V1: float, M1: float, T1: float, V2: float, M2: float, T2: float
    ):
        if V2 is not None and M2 is not None:
            T2: float = (M1 / M2) * T1 * (V2 / V1)
        elif V2 is not None and T2 is not None:
            M2: float = (V2 / V1) * M1 * (T1 / T2)
        elif T2 is not None and M2 is not None:
            V2: float = (M2 / M1) * V1 * (T2 / T1)
        else:
            raise ValueError(
                "For an isobaric process, out of Volume, Mass and Temperature, any 2 should be constrained"
            )

        return V2, M2, T2

    def __isochoric_final_state(
        self, P1: float, M1: float, T1: float, P2: float, M2: float, T2: float
    ):
        if P2 is not None and M2 is not None:
            T2: float = (M1 / M2) * T1 * (P2 / P1)
        elif P2 is not None and T2 is not None:
            M2: float = (P2 / P1) * M1 * (T1 / T2)
        elif T2 is not None and M2 is not None:
            P2: float = (M2 / M1) * P1 * (T2 / T1)
        else:
            raise ValueError(
                "For an isochoric process, out of Pressure, Mass and Temperature, any 2 should be constrained"
            )

        return P2, M2, T2

    def __isothermal_final_state(
        self, P1: float, M1: float, V1: float, P2: float, M2: float, V2: float
    ):
        if P2 is not None and V2 is not None:
            M2: float = (V2 / V1) * M1 * (P2 / P1)
        elif P2 is not None and M2 is not None:
            V2: float = (P1 / P2) * V1 * (M2 / M1)
        elif V2 is not None and M2 is not None:
            P2: float = (V1 / V2) * P1 * (M2 / M1)
        else:
            raise ValueError(
                "For an isothermal process, out of Pressure, Mass and Volume, any 2 should be constrained"
            )

        return P2, M2, V2

    def __adiabatic_final_state(
        self,
        P1: float,
        M1: float,
        V1: float,
        T1: float,
        P2: float,
        M2: float,
        V2: float,
        T2: float,
    ):
        y = self.__initial_state.gamma
        if P2 is not None and M2 is not None:
            V2 = ((P1 / P2) ** (1 / y)) * V1
            T2 = (M1 / M2) * T1 * (V1 / V2) ** (y - 1)
        elif P2 is not None and T2 is not None:
            V2 = ((P1 / P2) ** (1 / y)) * V1
            M2 = (T1 / T2) * M1 * ((V1 / V2) ** (y - 1))
        elif M2 is not None and T2 is not None:
            V2 = (M1 / M2) * (T1 / T2) * (V1 ** (y - 1))
            P2 = P1 * ((V1 / V2) ** y)
        elif M2 is not None and V2 is not None:
            P2 = ((V1 / V2) ** y) * P1
            T2 = (M1 / M2) * T1 * ((V1 / V2) ** (y - 1))
        else:
            raise ValueError(
                "For an adiabatic process, out of Pressure, Mass, Temperature and Volume, any 2 should be constrained"
            )

        return P2, V2, M2, T2

    @property
    def initial_state(self):
        return self.__initial_state

    @property
    def final_state(self):
        return self.__final_state
