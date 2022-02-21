# The transport library
ReactionEngine lets you calculate the viscosity, thermal conductivity and diffusion coefficients.
The code relies on *transport.dat* for the evaluation of these properties. A sample of the *transport.dat* file is shown below.

```
AR                 0   136.500     3.330     0.000     0.000     0.000
CH4                2   141.400     3.746     0.000     2.600    13.000
CO2                1   244.000     3.763     0.000     2.650     2.100
CO                 1    98.100     3.650     0.000     1.950     1.800
H2                 1    38.000     2.920     0.000     0.790   280.000
H2O                2   572.400     2.605     1.844     0.000     4.000
O2                 1   107.400     3.458     0.000     1.600     3.800
```

The different columns present in the above table is described below.

- Column-1: Name of species
- Column-2: Geometric configuration of the species; 0- single atom, 1- linear molecule, 2-nonlinear molecule
- Column-3: Lennard-Jones potential well depth expressed as $\epsilon/k_B$ and has units K
- Column-4: Lennard-Jones collision diameter ($\sigma$) in Angstroms
- Column-5: Dipole moment ($\mu$) in Debye
- Column-6: Polarizability ($\alpha$) in cubic Angstroms
- Column-7: Rotational relaxation collision number ($Z_\mathrm{rot}$) at 298 K


## Diffusion coefficients
### Binary diffusion coefficients
The binary diffusion coefficient is expressed as 
```math
D_{jk} = \frac{3}{16} \frac{ \sqrt{2\pi N_Ak_B^3T^3/m_{jk}} }{p\pi \sigma_{jk}^2 \Omega^{(1,1)} }
```
Here $N_A$ is the Avogadro's number, $k_B$ is the Boltzmann constant, $p$ is the pressure, and $T$ is the temperature. The reduced molar mass of the species pair (j,k) is defined as

```math
m_{jk} = \frac{m_jm_k}{m_j+m_k}
```

The collision integral $\Omega^{(1,1)}$ is determined using the reduced temperature $T^*_{jk}$ and the reduced dipole moment $\delta^*_{jk}$

```math
T^*_{jk} = \frac{k_BT}{\epsilon_{jk}}
```

```math
\delta^*_{jk} = \frac{1}{2} \frac{\mu^2_{jk}}{\epsilon_{jk}\sigma_{jk}^3}
```

The reduced dipole moment depends on the polarizability of the interacting molecules. In the case where both molecules are 
either polar or non-polar, then it follows

```math
\epsilon_{jk} = \sqrt{\epsilon_j \epsilon_k}
```

```math
\mu_{jk}^2 = \mu_j \mu_k
```

and the reduced collision diameter $\sigma_{jk}$ is defined as

```math
\sigma_{jk} = \frac{\sigma_k+\sigma_j}{2}
```

In the case of interaction between a polar and non-polar molecule

```math
\epsilon_{jk}= \xi^2  \sqrt{\epsilon_j \epsilon_k}
```

```math
\sigma_{jk} = \frac{1}{2} (\sigma_j + \sigma_k) \xi^{-1/6}
```

```math
\mu_{jk}= 0
```

```math
\xi = 1+\frac{1}{4} \alpha_n^* \mu_p^* \sqrt{ \frac{\epsilon_p}{\epsilon_n} }
```
Note that the subscripts $p$ and $n$ represents either $j$ or $k$. If $j$ is polar species, then the subscript $p$ refers to that species

```math
\alpha_n^* = \frac{\alpha_n}{\sigma_n^3}
```

```math
\mu_p^* = \frac{\mu_p}{\sqrt{\epsilon_p\sigma_p^3}}
```
The estimation of $\Omega^{(1,1)}$ is a table look up procedure using the values of $T^*_{jk}$ and $\delta^*_{jk}$.

### Mixture diffusion coefficients
The mixture diffusion coefficients are calculated from

```math
D_{k,m} = \frac{1-Y_k}{\sum_{j\ne k}^N X_j/D_{jk}}
```
In the above equation, $Y_k$ is the mass fraction of the species $k$, and $X_j$ is the mole fraction of species $j$. $N$ is the total number of species

## Viscosity
The viscosity of the mixture is calculated using pure species viscosity. The pure species viscosity is defined as

```math
\eta_k = \frac{5}{16} \frac{\sqrt{\pi m_k k_BT/N_A}}{\pi \sigma_k^2 \Omega^{(2,2)}}
```
The Lennard-Jones collision integral $\Omega^{(2,2)}$ is estimated using a table lookup procedure and depends on the reduced temperature. 

```math
T^*_k = \frac{k_BT}{\epsilon_k}
```
and the reduced dipole moment

```math
\delta^*_k = \frac{1}{2} \frac{\mu_k^2}{\epsilon_k \sigma_k^3}
```
The viscosity of the mixture is then defined as

```math
\eta = \sum_{k=1}^{N} \frac{X_k\eta_k}{\sum{j=1}^K X_j\Phi_{kj}}
```

```math
\Phi_{kj}= \frac{1}{\sqrt{8}} \left(1+\frac{M_k}{M_j}\right)^{-1/2} \left( 1+ \left(\frac{\eta_k}{\eta_j}\right)^{1/2} \left(\frac{M_j}{M_k}\right)^{1/4}\right)^2
```

## Thermal conductivity
Similar to calculating viscosity, the thermal conductivity of a mixture is calculated from the pure species thermal conductivity. 
The pure species thermal conductivity is defined as

```math
\lambda_k = \frac{\eta_k}{M_k} \left( f_t C_{v,t} + f_r C_{v,r} + f_v C_{v,v} \right)
```

```math
f_t = \frac{5}{2}\left( 1-\frac{2}{\pi} \frac{C_{v,r}}{C_{v,t}} \frac{A}{B} \right)
```

```math
f_v = \frac{\rho D_{kk}}{\eta_k}
```


```math
f_r = f_v \left( 1+ \frac{2}{\pi} \frac{A}{B}\right)
```

```math
A = 2.5-f_v
```

```math
B = Z_{rot} + \frac{2}{\pi} \left( \frac{5}{3} \frac{C_{v,r}}{R} + f_v\right)
```
The molar heat capacity $C_v$ for rotational, vibrational or transnational mode depends on the molecule's geometry. For linear molecules

```math
\frac{C_{v,t}}{R} = \frac{3}{2}
```
```math
\frac{C_{v,r}}{R} = 1
``` 
```math
\frac{C_{v,v}}{R} = C_v  - 2.5R
```
```math
C_v = C_p -R
```

For non-linear molecules
```math
\frac{C_{v,t}}{R} = \frac{3}{2}
```
```math
\frac{C_{v,r}}{R} = \frac{3}{2}
``` 
```math
\frac{C_{v,v}}{R} = C_v  - 3R
```
```math
C_v = C_p -R
```
