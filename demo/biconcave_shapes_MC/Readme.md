# Simplest Monte Carlo simulation
<!-- gToDo:   this doc needs to be updated-->
This folder contains the simplest Monte Carlo simulation that I could imagine, which still generates visual results that are non-trivial.

- [Theoretical outline of the simulation](#theoretical-outline-of-the-simulation)
- [Implementation](#implementation)
- [Data visualization](#data-visualization)

## Theoretical outline of the simulation
The goal of the simulation is to initiate a spherical triangulation and use a metropolis scheme to update the spherical triangulation while enforcing the following energy

![energy function](https://latex.codecogs.com/svg.latex?\Large&space;E=\frac{\kappa}{2}\int\mathrm{d}A(2H)^2+K_A\frac{(A-A_0)^2}{A_0}+K_V\frac{(V-V_t)^2}{V_t},)

where the first term is the Canham-Helfrich energy and minimizes the curvature (which makes the triangulation behave like a biological membrane [[Canham1970](https://doi.org/10.1016/S0022-5193(70)80032-7), [Helfrich1973](https://doi-org.tudelft.idm.oclc.org/10.1515/znc-1973-11-1209)], the second and third terms are fixing the area and volume respectively, with harmonic potentials. The area is fixed to its initial value, whereas the target volume is `60%` of the original sphere volume, i.e.

![target volume](https://latex.codecogs.com/svg.latex?\Large&space;V_t=0.6V_0.)

The constant ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;\kappa) is the bending rigidity determining how easy it is to bend the membrane. The constants ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_A) and ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_V) are Lagrange multipliers that fix the area and volume, respectively. Their values do not have physical interpretations, but the larger they are, the more important it becomes for the simulation to fix their respective values precisely. This means that if we make these constants too small, the deviations from desired target values will be unacceptably large. However, if we make these constants too large, the area and volume terms will penalize even small deviations too strongly, which means that all Monte Carlo steps will be rejected as soon as the target values are reached, and the simulation will ignore the curvature term.
Choosing ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_A) and ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_V) is usually done by trial and error (that's what I did here).

## Implementation

### Overview
There are four main ingredients that we need to implement with *flippy* to simulate the above-described system.
1. Implementation of the energy function and the `struct` containing all the function parameters.
2. Initiate a triangulation.
3. Initiate the built-in Monte Carlo updater and connect it with the initiated triangulation and the energy function.
4. Write the update loop. This step specifies in which order the nodes should be moved and flipped. Any external parameter should change during the updating.
   1. We will use this step to slowly transform the initial condition of the triangulation (a sphere) to a configuration that has the proper target volume,
   2. then, we will equilibrate the triangulation at that target volume.

### Energy function
The `MonteCarloUpdater` class of *flippy* will use the energy function we define here. This means that the signature of the function needs to be exactly what the `MonteCarloUpdater` class expects, namely:
```c++
double surface_energy(fp::Node<double, unsigned int> const&,
                      fp::Triangulation<double, int> const& ,
                      EnergyParameters const& )
```
The first argument of the function is the constant reference to the node (that is being updated), the second argument is a constant reference to the triangulation that the node is part of, and the first argument can be anything as long as we specify that in the `MonteCarloUpdater` declaration [later](#monte-carlo-updater-declaration). In this case, we want to pass a `struct` that holds all the parameters necessary for the energy function. We call this custom struct `EnergyParameters`, which is defined as follows:
```c++
struct EnergyParameters{double kappa, K_V, K_A, V_t, A_t;};
```
and later initialized as:
```c++
EnergyParameters prms{.kappa=10 /*kBT*/,
                          .K_V=100 /*kBT/area*/, .K_A=1000 /*kBT/volume*/,
                          .V_t=0.6*sphere_vol(R), .A_t=sphere_area(R)};
```
where
- `R` is the initial radius of the triangulation and is initialized such that the triangles barely not overlap
- `sphere_vol` and `sphere_area` are functions defined in the main file and return the volume and area of a sphere of a given radius.

Since this simulation does not have a physical length scale, the values for area and volume targets are given in arbitrary units.

The actual body of the energy function then looks as follows:
```c++
double surface_energy([[maybe_unused]]fp::Node<double, unsigned int> const& node,
                      fp::Triangulation<double, int> const& trg,
                      EnergyParameters const& prms){
    double V = trg.global_geometry().volume;
    double A = trg.global_geometry().area;
    double dV = V-prms.V_t;
    double dA = A-prms.A_t;
    double energy = prms.kappa*trg.global_geometry().unit_bending_energy +
                    prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return energy;
}
```
Where the input argument node has the compiler directive `[[maybe_unused]]` prepended, which tells the compiler that this variable is not used in the function body.

### Triangulation declaration
To create a spherical triangulation, we can initiate a `Triangulation` class:
```c++
fp::Triangulation tr(n_triang, R, r_Verlet);
```
here `fp` is *flippy*'s namespace. The first two template arguments specify the types of floating point and integral numbers the triangulation should use. One could have specified `float` and `short` instead. The third argument specifies the type of triangulation. In this case, we want to use a spherical one. However, the last template argument is optional since spherical triangulation is already the default. This means that we can also write:
```c++
fp::Triangulation guv(n_triang, R, r_Verlet);
```
The arguments of the triangulation instantiation are `n_triang`, which specifies the level of triangulation. *flippy*'s current implementation starts from an icosahedron, and after the `n_triang` steps of sub-triangulation, one gets

![node number](https://latex.codecogs.com/svg.latex?\Large&space;N_{nodes}=12+30n_{triang}+20\frac{n_{triang}(n_{triang}-1)}{2},)

number of nodes. The second variable, `R`, specifies the radius of the triangulated sphere, and the variable `r_Verlet` specifies the Verlet radius for the triangulation nodes. This sets the size of the [Verlet list](https://en.wikipedia.org/wiki/Verlet_list), which lists spatially close nodes. This information is necessary to implement the non-self-intersection property of the membrane efficiently.

After the above declaration, we will access the triangulation via the declared variable `guv`.

### Monte Carlo updater declaration
The Monte Carlo updater needs information about the triangulation, the energy function, and the parameter struct, all of which we have already defined and initiated. Additionally, we also need a random number generator.

The declaration can be made as follows:
```c++
fp::MonteCarloUpdater<double, int, EnergyParameters, std::mt19937,
                      fp::SPHERICAL_TRIANGULATION>
                      mc_updater(guv, prms, surface_energy, rng, l_min, l_max);
```
The template arguments again specify what type of updater we want.
- `double` and `unsigned int` specify the floating point and non-negative integral number types used in updating (same as in triangulation). We must use the same floating point type here as in the triangulation and the energy function's return value.
- `EnergyParameters` specifies the type of the third argument of the energy function.
- `std::mt19937` specifies the type of the random number generator. The provided value is a Mersenne Twister random number generator, part of the `c++` standard library.
- `fp::SPHERICAL_TRIANGULATION` specifies the triangulation type (not optional in this case), which has to match the triangulation type of the `Triangulation` class.

The instantiation parameters have the following meaning:
 - `guv` name of the `Triangulation` class instance we declared.
 - `prms` name of the `EnergyParameters` struct instance we declared.
 - `surface_energy` is the name of the energy function we defined.
 - `rng` is the name of the random number generator instance. This generator needs to be declared before the `mc_updater` in the code as
```c++
std::random_device random_number_generator_seed;
std::mt19937 rng(random_number_generator_seed());
```
 - `l_min` minimum distance between the triangulation nodes allowed during updating.
 - `l_max` maximum distance between connected triangulation nodes allowed during updating.

The instance `mc_updater` now provides access to functions that can attempt an update of the triangulation.

- `move_MC_updater` expects a node and a displacement vector and will attempt to update that node's position by the displacement vector
- `flip_MC_updater` expects a node and will randomly choose one of the neighbors of that node and attempt to flip the bond between the two nodes.


### The update loop
The simplest update loop would be one where we loop through all nodes and use the `MonteCarloUpdater` to attempt an update, then repeat this for a set number of times specified by `max_mc_steps`. This would look like this:

```c++
for(int t=1; t<max_mc_steps+1; ++t){
    for (auto const& node: guv.nodes()) {
      displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
      mc_updater.move_MC_updater(node, displ);
      mc_updater.flip_MC_updater(node);
  }
}
```
Here we used some new variables and functions that we defined before the loop in the main file:
- `displ` is a 3-dimensional vector of displacements of type `fp::vec3<double>`, which is a *flippy* builtin type.
- `displ_distr` is a uniform distribution from which displacements are drawn.

This update loop would work, but usually, we want to do a bit more than what's implemented there.
For starters, we want to separate the moving events from flipping events. We want to create a separate vector for node IDs which we can shuffle in each Monte Carlo step. And finally, we want to squish the initial sphere a bit to break the initial spherical symmetry since the final biconcave shapes that we expect are oblate. This we will do with the `scale_node_coordinates` method of triangulation.

```c++
tr.scale_node_coordinates(1, 1, 0.8); // squish the sphere in the z-direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly.

std::vector<int> shuffled_ids;
shuffled_ids.reserve(guv.size());
for(auto const& node: guv.nodes()){ shuffled_ids.push_back(node.id);}

for(int mc_step=0; mc_step<max_mc_steps; ++mc_step){
    for (int node_id: shuffled_ids) { // we first loop through all the beads and move them
        displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
        mc_updater.move_MC_updater(guv[node_id], displ); // guv[node_id] returns the node which has id=node_id
    }

    std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng); // then we shuffle the bead_ids

    for (int node_id: shuffled_ids) { // then we loop through all of them again and try to flip their bonds
        mc_updater.flip_MC_updater(guv[node_id]);
    }
}
```

### Data saving
The least effort way to save a snapshot of the triangulation is to use the built-in `make_egg` method. Which serializes every node of the triangulation to a `JSON` object.
```c++
fp::Json data_final = guv.make_egg_data();
```
Then we can use one of the *flippy*'s helper functions, `dump_json` to save the data to a file:
```c++
fp::json_dump("test_run_final", data_final);
```
The above command will save the data to a `test_run_final.json` file in the same folder where the executable was executed.

If one wants to continue the simulation after the final configuration, then the saved data can be used to initialize a triangulation:
```c++
fp::Json loaded_data = fp::json_read("test_run_final.json");
fp::Triangulation<double, int> loaded_guv(loaded_data, r_Verlet);
```
## Data visualization

This folder also provides a `data_vizualization.py` python file. If this file is run in the same folder, where
the `test_run_init.json` and `test_run_final.json` JSON files are saved (the executable output files), then the python file will generate two plots.
That of the initial and final configurations of a run.
