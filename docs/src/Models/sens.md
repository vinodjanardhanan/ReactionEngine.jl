# Sensitivity Analysis
The sensitivity analysis in ReactionEngine evaluates the global sensitivity using Sobol sampling. It reports the Sobol indices for all species mentioned in the input file against the micro-kinetic model parameters. The function evaluation is performed using a batch reactor model, and therefore, you may please refer to the documentation for the batch reactor model to understand the governing equations that are solved. 

In general, only the forward reaction rate parameters can be specified independently for a micro-kinetic model where the reactions are expressed as reversible pairs. The reverse reaction rate parameters are then calculated from the equilibrium constant. However, this requires thermodynamic data for the calculation of equilibrium constants. 
In many surface reaction mechanisms, the reversible reactions are explicitly expressed as two reactions; one represents the forward direction, and the other represents the reverse direction. This difficulty is mainly due to the non-availability of the NASA polynomial coefficients for the surface adsorbed species to calculate equilibrium constants. However, for such a reaction mechanism to be thermodynamically consistent, one has to ensure the following constraint on all forward-reverse reaction pairs.

```math
\frac{A_fT^\beta}{A_r T^\beta} = \exp(\frac{\Delta S}{R})
```
In such a case, all parameters correlate to at least one parameter. Therefore altering one parameter will alter another parameter constraint to the above correlation. While developing a surface reaction mechanism, one may assume that all parameters are uncorrelated. For instance, if the reaction mechanism contains 20 reversible reactions expressed as 40 irreversible pairs, one may want to perform a full sensitivity analysis assuming that all 40 parameters are uncorrelated. At some stage during the mechanism development, if you make the mechanism thermodynamically consistent, you may prefer to perform the analysis using only 20 independent parameters. In that case, altering any of the 20 parameters will also alter the parameter of the corresponding reaction pair. One can also think of secondary interactions, which are not considered here. So, the analysis assumes that the parameters are independent or uncorrelated in either case. This is important for performing Sobol decomposition or ANOVA decomposition. The Sobol decomposition of a function $y=f(x)$ is given as

```math
y= f(x) = f_0 + \sum_{i=1}^n f_i(x_i) + \sum_{i=1}^n \sum_{j=i+1}^n f_{ij}(x_i,x_j) + \ldots + f_{12\cdots n}(x_1, x_2,\ldots,x_n)
```
where
```math
f_0 = E(y)
```

```math
f_i(x_i) = E(y\vert x_i) - E(y)
```

```math
f_{ij}(x_i,x_j) = E(y\vert x_i,x_j) - E(y\vert x_i) - E(y\vert x_j) - E(y) = E(y\vert x_i,x_j) - f_i - f_j - f_0
```
Taking the variance of function $f(x)$

```math
V(y) = V[f_0] +V\left[ \sum_{i=1}^n f_i(x_i)\right] + V\left[\sum_{i=1}^n \sum_{j=i+1}^n f_{ij}(x_i,x_j)\right] + \ldots +V[ f_{12\cdots n}(x_1, x_2,\ldots,x_n)]
```
i.e.,

```math
V(y) = 0 +  \sum_{i=1}^n V[f_i(x_i)] + \sum_{i=1}^n \sum_{j=i+1}^n V[f_{ij}(x_i,x_j)] + \ldots + V[f_{12\cdots n}(x_1, x_2,\ldots,x_n)]
```

```math
V = \sum_{i=1}^n V_i + \sum_{i=1}^n \sum_{j=i+1}^n V_{ij}+ \ldots + V_{12\ldots n}
```

The global sensitivity index is then defined as

```math
S_i =\frac{ V_i}{V}
```
which is also known as the main sensitivity index of the first order Sobol index. The total sensitivity index is defined as

```math
S_i^T = \frac{V_i^T}{V}
```
where
```math
V_i^T= V_i + \sum_{j} V_{ij} + \sum_{jk} V_{ijk} + \ldots
```
The first order Sobol index gives the variance of the model output due to the corresponding parameter input. The total sensitivity index also considers the interaction among the parameters. Therefore, the difference between the total and main sensitivity indices shows the contribution due to interaction. 

## Running the code
The code is invoked by using the following method.
```@docs
global_sensitivity(input_file::AbstractString, lib_dir::AbstractString)
```
On the Julia REPL 
```julia
julia>using ReactionEngine
julia>global_sensitivity("sensitivity.xml","../lib/")
```

## Input file
The method takes *file\_path*  and path to *lib* directory as the arguments. The file_path points to the input XML file. The structure of the XML input file is shown below.

```
?xml version="1.0" encoding="ISO-8859-1"?>
<sensitivity>
    <gsa_species>CH4 H2O H2 CO CO2 O2 N2</gsa_species>
    <gsa_lb>5</gsa_lb>
    <gsa_ub>10</gsa_ub>
    <gsa_N>100</gsa_N>
    <gsa_model>cstr.xml</gsa_model>
    <p_smech>1:6=7:12,15=16,16:20</p_smech>
</sensitivity>
```
The meaning of different tags is listed below.

- <sensitivity>: Root xml tag for sensitivity analysis
- <gsa\_species>: The species that need to be monitored for the parameter sensitivity analysis. The code will produce the output of parameter sensitivity of only those species which are listed in this element
- <gsa\_lb>: lower percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the lower bound considered for design matrix generation is 90 with a gsa_lb = 10
- <gsa\_ub>: upper percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the upper bound considered for design matrix generation is 110 with gsa_ub = 10
- <gsa\_N>: the size of the sample size. With a sample size of 100, 100 points will be sampled between the lower and upper bounds for each parameter. The higher the sample size higher will be calculation time. 
- <gsa\_model>: the reactor model that will be considered for integration. The content of the tag is the name of the XML input file of the reactor model. The file must be present in the working directory. This is 
- <p\_smech>: parameters of the surface reaction mechanism that needs to be considered for sensitivity analysis. There are different options for specifying the parameter list from a mechanism
    -   All parameters: Use ":" as the value of <p\_smech> element to consider all pre-exponential factors and sticking coefficient as sensitivity parameters. In such specifications, the reverse reaction rates will not be adjusted to maintain the initial ratio. (e.g. <p\_smech>:</p_smech>)
    -   Reversible reaction: If the reverse reaction parameters need to be adjusted according to the initial ratio, then a reverse reaction must be specified for each forward reaction. For instance, if you want to consider the pre-exponential factor of reaction 15 as a sensitivity parameter, during sensitivity analysis, this parameter will be sampled between its upper and lower limit. This also changes the ratio between the forward and reverse reaction rate parameters as specified in the original mechanism. By specifying a reverse reaction, the ratio can be maintained a constant. For instance, in the above example, reaction 16 is the reverse of reaction 15. For every sample values of the pre-exponential factor of reaction 15, the pre-exponential factor of reaction 16 will be adjusted to maintain the initial ratio between these two in the mechanism input file. Only the parameter of reaction to the left of "=" will be considered as sensitivity parameter, and the pre-exponential factor to the right of "=" will be adjusted to maintain the initial ratio. 
    -   Contiguous parameters: Use ":" to specify contiguous reactions. For instance, 1:6=7:12 means pre-exponential factors of reactions from 1 to 6 will be considered as sensitivity parameters, and reactions from 7 to 12 are reverse pairs of reactions from 1 to 6. i.e., reaction 7 is the reverse pair of reaction 1, reaction 8 is the reverse pair of reaction 2 and so on. Comma-separated list of parameters can be provided in <p\_smech>
    -   Unconstrained reverse reactions: If the user does not want to constrain the reverse reaction parameters according to the initial ratio, then omit "=". For instance, <p\_smech>16:20</p_smech> means the pre-exponential factors of reactions 16 to 20 will be considered as sensitivity parameters. However, their corresponding reverse reaction parameters will not be adjusted to maintain the initial ratio.
