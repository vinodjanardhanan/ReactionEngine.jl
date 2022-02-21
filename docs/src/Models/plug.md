# Plug flow reactor
 The Plug flow reactor model simulates one dimensional tubular reactor. The governing equations solved are
 
 ```math
\frac{d (\rho u Y_k)}{dz} = \dot{s_k} M_k \frac{A_{s/L}}{A_c} + \dot{\omega_k}M_k 
```

```math
\frac{d (\rho u)}{dz} = \sum_{k=1}^{N_g} \dot{s_k} M_k \frac{A_{s/L}}{A_c} 
```

 In the above equations, $\rho$ is the density (kg/m$^3$), u is the velocity (m/s), $Y_k$ is the mass fraction of species $k$, z is the axial coordinate (m), $\dot{s_k}$ is the molar production rate of species $k$ due to surface reactions (mol/m$^2$-s), $\dot{\omega}_k$ is the molar production rate of species $k$ due to gasphase reactions (mol/m$^3$-s), $M_k$ is the molecular weight of species $k$, $A_{s/L}$ is the surface area per unit length (m), $A_c$ is the cross sectional area (m$^2$) and $\eta$ is surface area enhancement factor. This factor accounts for the actual surface area available for the surface reactions over the actual geometric surface area of the tubular reactor. 
 
 ## Running the code 
 
 
The code is invoked by using the following method.
```@docs
plug(input_file::AbstractString, lib_dir::AbstractString)
```
On the Julia REPL 
```julia
julia>using ReactionEngine
julia>plug("plug.xml","../lib/")
```
## Input file
The method takes *file\_path* as the argument. The file\_path points to the input XML file. The structure of the XML input file is shown below.

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<plug>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.25,CO2=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <length>0.3</length>
    <dia>0.005</dia>
    <u>0.1</u>
    <Tw>1073.15</Tw>
    <isothermal>true</isothermal>
    <cat-geom-factor>1000</cat-geom-factor>
    <surface_mech>data/ch4ni.xml</surface_mech>
</plug
```
The meaning of different tags is specified below.

- <plug> : The root XML tag for Plug
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab
- <molefractions> : inlet mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa
- <length> : length of the reactor in m
- <dia> : diameter of the reactor in m
- <u> : inlet velocity of the reacting mixture in m/s
- <Tw> : wall temperature in K. This option is provided for performing non-isothermal simulation, which is not supported in the current release
- <isothermal> : a boolean which accepts wither true or false. For the current release this must be true
- <cat-geom-factor> : surface area enhancement factor (refer $\eta$ in the governing equations)
- <surface_mech> : name of the surface reaction mechanism. Must be specified along with the path

# Output
The code generates two output files in the same directory where the input file **`cstr.xml`** is present. 
The file **`gas_profile.dat`** contains the mole fraction of the gas phase species as a function of time.
The file **`surf_profile.dat`** contains the surface coverages of adsorbed species as a function of time. 
In addition to these two files, the code also generates terminal output, which shows integration progress.
The terminal output is nothing by the integration time. 

An example terminal output is shown below
```
julia> plug("test/plug/plug.xml")
0.0000e+00
1.1015e-09
1.6078e-09
2.1141e-09
...
...
...
2.6191e-01
2.8298e-01
3.0000e-01
Solver integration: Success
```

A sample output of **`gas_profile.dat`** is shown below 
```
         z           T           p           u         rho         CH4         H2O          H2          CO         CO2          O2          N2
0.0000e+00  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  0.0000e+00  0.0000e+00  0.0000e+00  2.5000e-01  0.0000e+00  5.0000e-01
8.8802e-15  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  1.5726e-11  5.9064e-12  3.7358e-11  2.5000e-01  1.5402e-23  5.0000e-01
8.8811e-11  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  1.5763e-07  5.8187e-08  3.7344e-07  2.5000e-01  1.5149e-19  5.0000e-01
...
...
...
2.8298e-01  1.0732e+03  1.0000e+05  1.4643e-01  2.2212e-01  1.2200e-02  5.6965e-03  3.1137e-01  3.2276e-01  6.5032e-03  -6.7405e-19 3.4147e-01
3.0000e-01  1.0732e+03  1.0000e+05  1.4643e-01  2.2212e-01  1.2200e-02  5.6965e-03  3.1137e-01  3.2276e-01  6.5032e-03  -4.3287e-13 3.4147e-01
```

A sample output of **`surf_profile.dat`** is shown below 
```
         z           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)
0.0000e+00  1.0732e+03  8.9517e-01  2.1543e-03  1.0249e-01  1.7333e-09  2.2731e-08  3.4120e-06  1.6186e-04  6.7900e-06  7.3130e-06  1.8236e-17  5.9645e-11  5.3919e-10  3.8717e-10
8.8802e-15  1.0732e+03  8.9517e-01  2.1543e-03  1.0249e-01  1.7333e-09  2.2731e-08  3.4120e-06  1.6186e-04  6.7900e-06  7.3130e-06  1.8236e-17  5.9645e-11  5.3919e-10  3.8717e-10
...
...
...
2.8298e-01  1.0732e+03  4.1508e-01  1.5854e-01  5.2598e-05  3.9230e-11  6.8172e-06  5.8269e-07  4.2629e-01  5.5308e-07  2.7830e-05  6.5821e-11  1.5517e-11  9.7290e-11  4.8802e-11
3.0000e-01  1.0732e+03  4.1508e-01  1.5854e-01  5.2598e-05  3.9230e-11  6.8172e-06  5.8269e-07  4.2629e-01  5.5308e-07  2.7830e-05  6.5820e-11  1.5517e-11  9.7290e-11  4.8802e-11
```
