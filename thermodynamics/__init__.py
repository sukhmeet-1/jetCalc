"""
Modeling gaseous species,
their thermodynamic states and processes,
along with support for thermodynamic cycle analysis. This module includes

`phase`: Abstractions for synthesizing mixtures of gaseous species from the thermodynamic data (Source: NASA CEA).

`state`: Abstractions for defining thermodynamic states of the custom gaseous species built using `phase`. Built in support for the equation of state

`process`: Abstractions for transitioning from one thermodynamic state, built using `state` to another state. Support for various types of thermodynamic processes commonly encountered in analyses.

"""
