# Equilibrate
Equilibrate is a program for the calculation of equilibrium composition. There are different methods for calculating the equilibrium composition of a reacting mixture. For a system that consists of $N$ species and $K$ elements ($N > K$), the elements are conserved, and the classical approach involves the solution of $N+K$ non-linear equations. This approach calculates the moles of different species ($N$ numbers) at equilibrium and $K$ number of Lagrangian multipliers.  An alternate method solves $N-K$ non-linear equations involving equilibrium constants and $K$ equations describing component activities. The latter approach is adopted in the equilibrate program.

## Running the code

The code is invoked by using the following method.
```@docs
equilibrate(input_file::AbstractString, lib_dir::AbstractString)
```
On the Julia REPL 
```julia
julia>using ReactionEngine
julia>equulibrate("equil.xml","../lib/")
```


## Input file
The method takes *file\_path* as the argument. The file_path points to the input XML file. The structure of the XML input file is shown below.

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<equil>
    <gasphase>CH4 H2 CO CO2 H2O O2 N2 </gasphase>
    <molefractions>CH4=0.6, H2O=0.1, CO2=0.2, N2=0.1</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
</equil>
```
The meaning of different tags is specified below.

- <equil> : The root XML tag for equilibrate
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab
- <molefractions> : Initial mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa

## Output
The code generates only screen output. An example output is shown below.

```
Inititial condition:

Species          moles           molefraction
       CH4       6.7244e+00      6.0000e-01 
       CO2       2.2415e+00      2.0000e-01 
       H2O       1.1207e+00      1.0000e-01 
        N2       1.1207e+00      1.0000e-01 
        CO       0.0000e+00      0.0000e+00 
        O2       0.0000e+00      0.0000e+00 
        H2       0.0000e+00      0.0000e+00 

Equilibrium composition @ T= 1073.15 K and p=100000.0 Pa

Species          moles           molefraction
       CH4       3.3912e+00      1.8973e-01
       CO2       1.2649e-02      7.0766e-04
       H2O       1.6324e-02      9.1329e-04
        CO       5.5621e+00      3.1118e-01
        O2       3.2984e-23      1.8454e-24
        H2       7.7709e+00      4.3476e-01
        N2       1.1207e+00      6.2703e-02
```        
