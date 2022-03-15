# Simplest Monte Carlo simulation
This folder contains the simplest monte carlo simulation that I could imagine building, which still generates visual results, that are non-trivial. 

- [Theoretical outline of the simulation](#theoretical outline of the simulation)
- [Implementation](#implementation)

## Theoretical outline of the simulation
The goal of the simulation is to initiate a spherical triangulation and use a metropolis scheme to update the spherical triangulation, while enforcing the following energy

![energy function](https://latex.codecogs.com/svg.latex?\Large&space;E=\frac{\kappa}{2}\int\mathrm{d}A(2H)^2+K_A\frac{(A-A_0)^2}{A_0}+K_V\frac{(V-V_t)^2}{V_t},)

Where the first term is the Canham-Helfrich energy and minimizes the curvature (which makes the triangulation behave like a biological membrane [[Canham1970](https://doi.org/10.1016/S0022-5193(70)80032-7), [Helfrich1973](https://doi-org.tudelft.idm.oclc.org/10.1515/znc-1973-11-1209)], the second and third terms are fixing the area and volume respectively, with harmonic potentials. The area is fixed to its initial value whereas the target volume is `60%` of the original sphere volume, i.e.

![target volume](https://latex.codecogs.com/svg.latex?\Large&space;V_t=0.6V_0.)

The constant ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;\kappa) is the bending rigidity which determines how easy it is to bend the membrane. The constants ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_A) and ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_V) are lagrange multipliers that fix the area and volume respectively. Their values do not have physical interpretations but the larger they are more important it becomes for the simulation to fix their respective values precisely. This means that if we make these constants too small the deviations from desired target values will be unacceptably large. However, if we make these constants too large, then the area and volume terms will penalize even small deviations too strongly, which means that all Monte Carlo steps will be rejected as soon as the target values are reached and the simulation will ignore the curvature term. 
Choosing ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_A) and ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_V) is usually done by trial and error (that's what I did here).

# Implementaion