
# Thermodynamic properties
The thermodynamic properties are calculated using the NASA polynomials (7 coefficient fit).
The library module that reads the NASA polynomials expect the input file "*therm.dat*" in the "**data**" folder
by default. The location for the "**data**" folder is at the same level as that of the "*src*" folder. The
NASA polynomial file can be downloaded from a number of websites. The library lets you calculate the
thermochemical data of pure species or for a gas mixture. A skeleton of the NASA polynomial is given below
```
 H2               TPIS78H   2   00   00   00G   200.000  3500.000   1000.00    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01 8.46810200E+03    4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
```


## Pure species propoerties
The pure species thermodynamic properties are
calculated according to the following equations. 

### Specific heat
The molar specific heat capacity $C_{pk}^0$ (J/mol-K) of a chemical species k is calculated using 
```math
\frac{C_{pk}^0}{R} = a_{1k} + a_{2k} T_k + a_{3k}T^2 + a_{4k}T^3 + a_{5k} T^4
```

```@docs
IdealGas.cp
```

### Enthalpy
The molar enthalpy $H_{pk}^0$ (J/mol) of a chemical species k is calculated according to 
```math
\frac{H_k^0}{RT} = a_{1k} + \frac{a_{2k}}{2} T + \frac{a_{3k}}{3}T^2 + \frac{a_{4k}}{4}T^3 + \frac{a_{5k}}{5} T^4 + \frac{a_{6k}}{T}
```

```@docs
IdealGas.H
```

### Entropy
The molar entropy $S_{pk}^0$ (J/mol-K) of a chemical species kis calculated according to 
```math
\frac{S_k^0}{R} = a_{1k} \ln T + a_{2k} T + \frac{a_{3k}}{2}T^2 + \frac{a_{4k}}{3}T^3 + \frac{a_{5k}}{4} T^4 + a_{7k}
```
```@docs
IdealGas.S
```

## Mixture average properties
The mixture average specific heat, enthalpy and entropy may be calculated using the following equations

### Specific heat
```math
\overline{C}_p = \sum_{k=1}^K C_{pk}X_k
```
```@docs
IdealGas.cpmix
```

### Enthalpy
```math
\overline{H} = \sum_{k=1}^K H_kX_k
```
```docs
IdealGas.hmix
```
### Entropy
```math
\overline{S} = \sum_{k=1}^K \bigg(S_k^0 -R\ln X_k - R\ln(p/P_\mathrm{atm}) \bigg)X_k.
```

### Gibb's free energy
```math
\overline{G} = \sum_{k=1}^K \bigg[ H_k - T_k\bigg(S_k^0 -R\ln X_k - R\ln(p/P_\mathrm{atm}) \bigg) \bigg]X_k
```