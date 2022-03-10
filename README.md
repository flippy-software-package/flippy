[![pipeline status](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/pipeline.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![coverage report](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/coverage.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases)
[![licence](https://img.shields.io/badge/licence-MIT-green)](https://gitlab.tudelft.nl/idema-group/flippy/-/blob/master/LICENSE)
[![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org)
# flippy

<img src="assets/flippy.png" alt="flippy" width="300"/> 

c++ package for dynamically triangulated membrane simulations.

[[_TOC_]]

# Gallery

<img src="https://surfdrive.surf.nl/files/index.php/s/6HCtX7B4NwE6w48/download" alt="rbc" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/1O3oJBNMT0vLxIv/download" alt="bubble_collision" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/mua39VAvdxmobD4/download" alt="cell division" width="300">

# Support
Flippy is still in active development and the documentation is almost non-existent, but I try my best to reply to e-mails and provide support for scientists who want to use flippy.
### for questions about general usage
please use the support email [![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org).
### for bugfixes 
please create an [issue](https://gitlab.tudelft.nl/idema-group/flippy/-/issues).
### for feature requests
again the [issues](https://gitlab.tudelft.nl/idema-group/flippy/-/issues) page can be used but be aware that new features will be slow to come.
### want to join the team?
write me an email!


# How to get it

*flippy* is a headers only library, so all you need to is to download the `flippy` subfolder and copy it into your project.

Or if you prefer using a single header file you can just download the [flippy.hpp](https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/single_header_flippy/flippy.hpp) header from the `single_header_flippy` folder.

# Documentation
You can find flippy's [user manual](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/User_Manual) and automatically generated code [documentation](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/Documentation) over on the [wiki](https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/docs/mainpage.md).

# Examples of usage

This is a simple example of a monte carlo update of a spherical triangulation in under 100 lines. The code will run for several seconds to several minutes depending on the strength of your cpu and will produce a biconcave shape, i.e. something resembling a red blood cell.
This code assumes that the single header `flippy.hpp` is in the same folder as the main file.
This code saves two data files in the folder where the binary file will be executed, under the names `test_run_init.json` and `test_run_final.json` which are full snapshots of the initial and final configurations.

This code can be found in the subfolder `demo/simplest_MC` together with a python script that visualizes the data, and a simple `CMake` file that can be used to build the code.
```cpp
//demo/simplest_MC/main.cpp

#include <random> // needed for random displacement generation
#include <vector> // need for std::vector
#include "flippy.hpp"

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easily passed to functions
    double bending_rigidity, K_V, K_A, R, V_t, A_t, linear_displ, l_max_square;
    int n_trng, max_mc_steps, surface_updates_per_mc_step;
}tension;

// This is the energy function that is used by flippy's built in updater to decide if a move was energetically favourable or not
double surface_energy_area_volume_ensemble([[maybe_unused]]fp::Node<double, int> const& node, fp::Triangulation<double, int> const& triangulation, SimulationParameters const& prms){
    double V = triangulation.global_geometry().volume;
    double A = triangulation.global_geometry().area;
    double dV = V-prms.V_t;
    double dA = A-prms.A_t;
    double e_tot = prms.bending_rigidity*triangulation.global_geometry().unit_bending_energy + prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return e_tot;
}

int main(){
    int n_triang = 7; // triangulation iteration number of nodes N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
    double l_min = 2;
    double R = l_min/(2*sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.))); // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length are close to minimal. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 2.*l_min; // if you make l_max closer to l_min bond_flip acceptance rate will go down
    double V0 = sphere_vol(R);
    SimulationParameters prms{
        .bending_rigidity = 10 /*kBT_*/, .K_V = 100 /*kBT_*/, .K_A=1000, .R=R /*a.u.*/,
        .V_t=V0, .A_t=sphere_area(R),
        .linear_displ=l_min/8., // side length of a voxel from which the displacement of the node is drawn
        .n_trng=n_triang,
        .max_mc_steps=1000, // max number of iteration steps (depending on the strength of your cpu, this should take anywhere from a couple of seconds to a couple of minutes
        .surface_updates_per_mc_step=50
    };

    std::random_device random_number_generator_seed;
    std::mt19937 rng(random_number_generator_seed()); // create a random number generator and seed it with current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, int> tr(prms.n_trng, prms.R, 2*l_min);
    fp::MonteCarloUpdater<double, int, SimulationParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION> mc_updater(tr, prms, surface_energy_area_volume_ensemble, rng, l_min, l_max);

    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector
    std::uniform_real_distribution<double> displ_distr(-prms.linear_displ, prms.linear_displ); //define a distribution from which the small displacements in x y and z directions will be drawn

    tr.scale_node_coordinates(1,1,0.8); // squish the sphere in z direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly
    fp::Json data_init = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_init", data_init);  // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<int> shuffled_ids;
    shuffled_ids.reserve(tr.size());
    for(auto const& node: tr.nodes()){ shuffled_ids.push_back(node.id);} //create a vector that contains all node ids. We can shuffle this vector in each MC step, to iterate randomly through the nodes
    double progress=0;
    double volume_adaptation_fraction = 0.3; // this fraction of the simulation progress is devoted to going from initial volume to target volume
    double cooldown_fraction = 0.1; // this fraction of the simulation progress is devoted to temperature cooldown

    for(int t=1; t<prms.max_mc_steps+1; ++t){
        if(t%10==0){
            std::cout<<"V/V_t="<<tr.global_geometry().volume/prms.V_t<<" A/A_t="<<tr.global_geometry().area/prms.A_t<<'\n';
            tr.make_global_geometry(); // Every so often we need to recalculate the global geometry to reduce the error accumulation, which happens because we constantly add to and subtract from the global geometry
            tr.make_verlet_list(); // We need to update the Verlet list from time to time since the local neighbourhood of the nodes changes
        }
        progress = t/((double)prms.max_mc_steps);
        if(progress<volume_adaptation_fraction){
            // for the first half of the simulation we gradually decrease the target volume. If we do not do this the volume term will dominate the energy, and we will get weird shapes!
            // since the MonteCarloUpdater takes a reference to the SimulationParameters struct, we can just change the parameters and flippy will have access to the updated version.
            prms.V_t = V0*(1.-0.4*(progress/volume_adaptation_fraction));
        }
        else if((progress>1.-cooldown_fraction)){
            mc_updater.reset_kBT(1. - (progress - volume_adaptation_fraction)); //here we cool down the temperature to force the monte carlo simulation to go to the minimum configuration faster
        }

        for(int update_step=0; update_step<prms.surface_updates_per_mc_step; ++update_step) {
            for (int node_id: shuffled_ids) { // we first loop through all the beads and move them
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                mc_updater.move_MC_updater(tr[node_id], displ);
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng); //then we shuffle the bead_ids
            for (int node_id: shuffled_ids) { // then we loop through all of them again and try to flip their bonds
                mc_updater.flip_MC_updater(tr[node_id]);
            }
        }

    }

    // MonteCarloUpdater counts the number of accepted and rejected moves, and it distinguishes between if a rejectio occurred because of the energy or the bond length constraint.
    // We can use this to print a simple statistics here. This will help us decide if our displacement size is too large for example.
    std::cout<<"percentage of failed moves: "<<(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((long double)mc_updater.move_attempt_count())<<'\n';
    std::cout<<"percentage of failed flips: "<<(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((long double)mc_updater.flip_attempt_count())<<'\n';

    fp::Json data_final = {{"nodes", tr.make_egg_data()}}; // tr.make_egg_data() is sufficient to create a json object, but the python file expects a nested structure with nudes key containing the state of the triangulation`
    fp::json_dump("test_run_final", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable
    return 0;
}
```

# Versioning

Crrent release [![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases) is experimental. This means that the API may change significantly. Every input on the usability of `flippy`'s API is very welcome.

This repository follows [Semantic Versioning](https://semver.org/) guidelines.

Given a version number *MAJOR*.*MINOR*.*PATCH*, the flippy project will increment the:

- *MAJOR* version when we make incompatible API changes,
- *MINOR* version when we add functionality in a backwards compatible manner, and
- *PATCH* version when we make backwards compatible bug fixes.

Additional labels for pre-release and build metadata are available as extensions to the *MAJOR*.*MINOR*.*PATCH* format.

as long as *MAJOR* version is 0 the API is unstable and any *MINOR* update can be backwards incompatible!

## well tested part of the api

- spherical triangulation
- vec3
- nodes
- debug utils / utils
- MonteCarloUpdater

## new and poorly tested

- planar triangulation

## coming soon

- tubular triangulation

## could be implemented at some point

- solid bodies
    - other objects that the triangulations could interact with
- force based updater
    - a utility class like MonteCarloUpdater, which uses force balance functions to update node positions

## changes in [![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases)

### breaking changes 
  - renamed `global_geometry.dA_K2` and `node.scaled_curvature_energy` to `unit_bending_energy`
    - `unit_bending_energy` also differs from `scaled_curvature_energy` by a factor of `0.5`
  - removed `debug_utils` from flippy. This functionality was unrelated to membrane simulations and simply offered additional printing and timing capabilities.
### new features
- none
### bugfixes
  - removed default constructor from `MonteCarloUpdater` since it was implicitly deleted anyway.
  - changed update counters types in `MonteCarloUpdater` to long instead of Index to avoid integer overflow.
  - double check if the stdlib defines M_PI and define it if not.

# licence

*flippy* is under MIT License, which means you can do pretty much anything with it. For more information read the `LICENCE` file.
