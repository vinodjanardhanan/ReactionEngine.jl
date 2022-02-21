
```@meta
CurrentModule = ReactionEngine
```
Documentation for [ReactionEngine](https://github.com/vinodjanardhanan/ReactionEngine.jl).

# ReactionEngine
Reaction Engine is a Julia package for the simulation of chemical reactors where gas-phase reactions and/or surface reactions occur. The package provides i) a library for calculating reaction source terms, transport properties, thermodynamic properties, and ii) a bunch of  reactor models for reactor simulation. Current capabilities are limited to the simulation of gas-phase and surface chemistry. A minimal description of equations used for the calculation of properties and the governing equations used for the simulation of different reactor models are provided in this manual. The current version has the following reactor models

- inspect: A program for integrating surface coverages
- thermo_probe : A program for calculating the thermochemical properties 
- transport_properties : A program for the calculation of transport properties (viscosity, thermal conducitivity, and diffusion coefficients)
- batch : A program for the simulation of batch reactor
- cstr : A program for the simulation of continuous stirred tank reactor
- plug : A program for the simulation of tubular flow reactor 
- gsa : A program for evaluating the global sensitivity coeffieint of reaction rate parameters
- equilibrate : A program for calculating equilibrium compositions