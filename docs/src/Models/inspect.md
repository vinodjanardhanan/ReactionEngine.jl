# Inspect
Inspect is a module for testing a surface reaction mechanism. The code integrates the surface coverages until a steady state is reached for a given gas phase composition, temperature and pressure. i.e.,

```math
\frac{d\theta_k}{dt} = \frac{\dot{s}_k \sigma_k}{\Gamma}
```
where $\theta_k$ is the surface coverage of species  $k$, $t$ is the time, $\dot{s}_k$ is the molar production rate (mol/m$^2$-s) of species $k$, $\sigma_k$ is the site coordination number and $\Gamma$ is the total site density (mol/m$^2$). After successful integration, the code returns the steady-state coverages, surface flux (reaction rate) of individual surface species and the sum of surface coverages and fluxes. Ideally, the surface coverages must sum to 1, and the sum of fluxes must vanish at steady state.

## Running the code

The code is invoked by using the following method.
```@docs
inspect(input_file::AbstractString, lib_dir::AbstractString)
```
On the Julia REPL 
```julia
julia>using ReactionEngine
julia>inspect("inspect.xml","../lib/")
```
## Input file
The method takes *input\_path* and path to the *lib* directory as the argument. The file_path points to the input XML file. The structure of the XML input file is shown below.

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<inspect>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.2,H2O=0.1,CO2=0.7</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <mech>data/ch4ni.xml</mech>
</inspect>
```
The XML element <gasphase> lists all gas-phase species participating in the surface reaction mechanism. The surface reaction mechanism is an XML file specified using the XML element <mech>. Please refer to the Library section for more information on mechanism input files. In the above example, the name of the surface reaction mechanism file is *ch4ni.xml* and is present in the data directory. The XML input for the mechanism may be present in any directory and not necessarily in the data folder. The *data* folder may be found at the same level as the *src* folder in the current distribution. It is advisable to maintain all library-related files in the data folder. 

In addition to the XML input file, the code also requires a *therm.dat* file. The code expects this file in the *data* folder by default. You must ensure that all the species present in the XML element <gasphase> are listed in the *therm.dat* file.

## Output
The code does not generate any file output.  An example of terminal output that inspect generates is shown below

```
Initial coverage: 
-----------------------
        (NI)     6.0000e-01 
       H(NI)     0.0000e+00 
       O(NI)     0.0000e+00 
     CH4(NI)     0.0000e+00 
     H2O(NI)     4.0000e-01 
     CO2(NI)     0.0000e+00 
      CO(NI)     0.0000e+00 
      OH(NI)     0.0000e+00 
       C(NI)     0.0000e+00 
     HCO(NI)     0.0000e+00 
      CH(NI)     0.0000e+00 
     CH3(NI)     0.0000e+00 
     CH2(NI)     0.0000e+00 

T(K):   1073.15
p(Pa):  100000.0

Final coverage and flux after [10.0](s): 

     Species     coverage        flux(mol/m2-s)
-------------------------------------------------
        (NI)     2.9259e-01      -2.2798e-14 
       H(NI)     3.4172e-03      +2.0295e-14 
       O(NI)     7.0353e-01      +4.2465e-14 
     CH4(NI)     4.5279e-10      -2.3961e-17 
     H2O(NI)     8.4369e-05      +7.3585e-16 
     CO2(NI)     8.5809e-06      +0.0000e+00 
      CO(NI)     1.4316e-04      +1.3353e-19 
      OH(NI)     2.2621e-04      -1.9896e-14 
       C(NI)     1.1795e-06      -4.2775e-19 
     HCO(NI)     3.0945e-16      -3.3087e-24 
      CH(NI)     1.8607e-10      +6.9880e-19 
     CH3(NI)     1.7535e-09      +7.9791e-19 
     CH2(NI)     1.2266e-09      -9.2538e-19 

Sum of fluxes (mol/m2-s): 2.0778686468957624e-14
Sum of coverages: 0.9999999989545821
```

