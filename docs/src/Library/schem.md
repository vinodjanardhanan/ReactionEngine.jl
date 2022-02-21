# Surface chemistry
## Input file
The input file for the surface chemistry is an XML file with the following structure, and the file is expected in the “**data**” folder at the same level as "*src*" folder. The outermost tag is *<surface_chemistry>*, which takes in two arguments. The argument *unit* takes in the unit of activation energy, and the attribute *name* accepts a string, which is the name of the mechanism. One of the limitations of the current implementation is that it expects the stoichiometric coefficient of participating species to be one. 

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<surface_chemisrty unit="kJ/mol" name="ch4ni">
    <species>(ni) H(ni) O(ni) CH4(ni) H2O(ni) CO2(ni) CO(ni) OH(ni) C(ni) HCO(ni) CH(ni)  CH3(ni) CH2(ni) </species>
    <site name="(ni)">
        <coordination>ch4(ni)=1,co(ni)=2.0</coordination>
        <density unit="mol/cm2">2.66e-09</density>
        <initial>h2o(ni)=0.4,(ni)=0.6 </initial>
    </site>
    <stick>
        <rxn id="1" >h2 + (ni) + (ni) => h(ni) + h(ni)   @ 1.0000e-2 </rxn>
        <rxn id="2" >o2 + (ni) + (ni) => o(ni) + o(ni)   @ 1.0000e-2 </rxn>
        ...
        ...
        <rxn id="6" >co + (ni) => co(ni)                 @ 5.0000e-1 </rxn>
    </stick>
    <arrhenius>
        <rxn id="7"  >h(ni) + h(ni) => (ni) + (ni) + h2  @ 2.545e+19  0.0  81.21  </rxn>
        <rxn id="8"  >o(ni) + o(ni) => (ni) + (ni) + o2  @ 4.283e+23  0.0  474.95 </rxn>
        ...
        ...
        <rxn id="42" >c(ni) + oh(ni) => o(ni) + ch(ni)   @ 1.625e+21  0.0  128.61 </rxn>
    </arrhenius>
    <coverage id="12 20 21">co(ni)=-50</coverage>
    <coverage id="23">co(ni)=50</coverage>
    <mwc>3 4</mwc>
    <order id="23 ">co(ni)=2</order>
    <order id="12 ">c(ni)=1.5</order>
</surface_chemisrty>
```
The order in which the tags appear in the XML file is immaterial. However, the main tags can appear only once. 

### Surface adsorbed species 

The XML element *<species>* specifies all the surface species, including the active site that participates in the reactions. The user may specify the names in upper case or lower case.  

### Site specific information
The XML tag *<site>* is used to specify the site-specific information such as active site, the site coordination number for the adsorbed species, the surface site density, and the initial coverages. An example snippet is shown below.

```
<site name="(ni)">
    <coordination>ch4(ni)=1,co(ni)=2.0</coordination>
    <density unit="mol/cm2">2.66e-09</density>
    <initial>h2o(ni)=0.4,(ni)=0.6 </initial>
</site>
```
The *<site>* tag has an attribute *name*, which is the name of the active site, and must be part of the *<species>* list. The tag *<coorination>* is used to specify the site coordination number for each of the adsobed species. By default, the site coordination number for all adsorbed species is one. The user need not specify the *<coordination>* tag if the site coordination number for all adsorbed species is one. However, if any adsorbed species occupies more than one adsorbed site, then those species must be specified in the *<coordination>* tag with their corresponding site coordination number. The site density is specified using the XML element *<density>*. The unit for site density can be either mol/m$^2$ or mol/cm$^2$. The initial state of the surface is specified using the XML element *<initial>*. The sum of coverages mentioned in the *<initial>* elements must add up to 1.0. The specification of this initial coverage has no bearing on the calculation results. The state of the surface is decided by the gas phase and the prevailing temperature and pressure conditions. The values specified in the *<initial>* element are only for initialization. 

### Sticking reactions

The adsorption reactions are specified as sticking reactions. All such reactions are specified within the XML tag *<stick>* with the XML element *<rxn>*. The *<rxn>* element takes *id* as an attribute, which must be unique. Ideally, the *id* must be a series representing the number of reactions present in the mechanism. The content of the *<rxn>* element is the reaction’s stoichiometry, followed by the reaction parameters. The stoichiometry and the parameters are separated by the ‘@‘ symbol. A ‘=>’ symbol in the stoichiometry represents a forward reaction. If reverse reactions need to be specified, the user must ensure that the NASA polynomials are available for the surface adsorbed species. Three parameters can be specified after the ‘@‘ symbol. The first parameter specifies the sticking coefficient, the second is the temperature dependency of the sticking coefficient, and the third is the activation energy. Generally, sticking reactions are specified only the sticking coefficient, and in case the other two parameters are specified, then the rate constant will be modified accordingly. For more details on how the sticking coefficient is converted to rate constant, please refer to the section on rate calculation. An example of the specification of the sticking reaction is shown below. 
```
<stick>
    <rxn id="1"> h2 + (ni) + (ni) => h(ni) + h(ni) @ 1.0000e-2 </rxn>
</stick>
```
In the above example, 1.0e-2 represents the sticking coefficient $s^0$. 

### Arrhenius reaction

The desorption reactions, dissociation reactions, and disproportionation reactions are specified in the Arrhenius form. These reactions are defined in the XML tag *<arrhenius>*. The *<rxn>* element is the same as that of the sticking reaction; however, the first reaction parameter, in this case, represents the pre-exponential factor for the reaction. The second and third represent the pre-exponential factor’s temperature dependence and the activation energy. The unit of the activation energy is specified as an attribute for the outermost element *<surface_chemistry>*. An example of Arrhenius reaction specification is shown below. 
```
<arrhenius>
        <rxn id="7"> h(ni) + h(ni) => (ni) + (ni) + h2  @ 2.545e+19     0.0     81.21   </rxn>
</arrhenius>
```
2.545e+19 is the pre-exponential factor in the above example, 0.0 is the temperature exponent $\beta$, and  81.21 is the activation energy in kJ/mol. The units of the pre-exponential factors are assumed as **cm-mol-s**. It does not support any other units currently. 

## Coverage dependent reactions
The activation energy of some of the reactions may depend on the surface coverage of adsorbed species. 
The coverage dependency of a reaction is expressed using the XML element *<coverage>*. The element takes an attribute *id*, a single string of reaction id’s separated by white space. For example, the following specification means,
```
<coverage id="12 20 21">co(ni)=-50</coverage>
```
Reactions 12, 20, and 21 are dependent on the surface coverage of CO. The modified activation energy, in this such cases will be

```math
E_i = E_i + \epsilon_{ki} \theta_{k},
```
Where $E_i$ is the activation energy of i'th reaction in the mechanism, $\epsilon_{ki}$ is the coverage dependent activation energy of i'th reaction on k'th species, and $\theta_k$ is the surface coverage of species k. All reactions that are dependent on the coverage of a particular adsorbed species must be specified using a single XML element *<coverage>*. A reaction can have multiple coverage dependencies. In such cases, multiple *<coverage>* elements must be used to specify the coverage dependency. For e.g. if reaction 2 is dependent on the coverage of CO by -50 kJ/mol and O by 20 kJ/mol, then its coverage dependency is represented as

```
<coverage id="2">co(ni)=-50</coverage>
<coverage id="2">o(ni)=20</coverage>
```

## Modification of reaction order
Since the mechanism assumes unit stoichiometric coefficients for the participating reactants and products, the order of the reaction with respect to
any reacting species will be one. Additional order dependency may be introduced using the XML element *<order>*. An example is shown below.

```
<order id="23 ">co(ni)=2</order>
```
This means in reaction 23, the order of the reaction w.r.t surface coverage of CO is 2. For such reactions, the rate constant will be multiplied by $\theta_k^{\mu_{ki}}$. Similar to the case of coverage dependency, multiple order dependence may be specified for a given reaction. An example is given below

```
<order id="23 ">co(ni)=2</order>
<order id="23 ">o(ni)=3</order>
```


Please refer to the rate calculation section for more details. 

## Calculation of reaction rates

The reaction rates are calculated using mean-field approximation. The relevant equations are given below.

### Conversion of sticking reactions

The sticking reactions are converted to standard rate constants according to 

```math
k = \frac{s^0_k}{\Gamma^\mathrm{n}} \sqrt{\frac{RT}{2\pi M_k}}
```
n is the number of free sites available for collision. $\Gamma$ is the total site density, $s^0_k$ is the sticking coefficient of species k, and $M_k$ is the molecular weight of the species. The sticking coefficient 
may also be expressed in terms of arrhenius expression according to

```math
s^0_k = aT^b \exp(-c/RT)
```


If Motz-Wise correction is specified for the reaction, then k is modified according to 

```math
k = \frac{k}{(1-0.5s^0_k)}
```

If additional parameters such as the temperature dependence of the pre-exponential factor and activation energy are specified for a reaction, then the rate constant becomes

```math
k = k T^\beta \exp(-E/RT)
```

The forward reaction rate constant for the $i$'th reaction in the mechanism is calculated according to 
```math
k_{fi} = k^0_i T^{\beta_i} \exp\left(-\frac{E_i}{RT}\right) \prod_{k=N_g+1}^{N_g+N_s} \theta_k^{\mu_{ki}} \exp\left( -\frac{\epsilon_{ki}\theta_k}{RT}\right)
```
Here $k^0_i$ is the pre-exponential factor for the reaction i $\beta_i$ is the power for the temperature dependence of the pre-exponential factor, $N_g$ and $N_s$ are respectively the number of gasphase and surface species. $\theta_k^{\mu_{ki}}$ is the term accounting for the modification of the reaction order as explained earlier, $\epsilon_{ki}$ is the coverage dependent activation energy of reaction i on species k, and $\theta_k$ is the surface coverage of species k. The rate of production or consumption of a species k is then calculated according to 
```math
\dot{s}_k = \sum_{i=1}^{N_r} \nu_{ki}  \left(  \prod_{k=1}^{N_g + N_s} k_{fi} [X_k]^{\nu_{ki}^{\prime}} - \prod_{k=1}^{N_g + N_s} k_{ri} [X_k]^{\nu_{ki}^{\prime\prime}}\right) 
```
where $\nu_{ki}^{\prime}$ is the stoichiometric coefficient of the species k on the reactant side and $\nu_{ki}^{\prime\prime}$ is the stoichiometric coefficient of species k on the product side, and $[X_k]$ is the surface concentration of species k. The relationship between the surface coverages and surface concentrations is given by
```math
\theta_k = \frac{[X_k]\sigma_k}{\Gamma}
```
where $\sigma_k$ is the site coordination number for species k.

