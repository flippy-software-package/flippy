[![pipeline status](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/pipeline.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![coverage report](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/coverage.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases)
[![licence](https://img.shields.io/badge/licence-MIT-green)](https://gitlab.tudelft.nl/idema-group/flippy/-/blob/master/LICENSE)
[![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org)
# flippy
c++ package for dynamically triangulated membrane simulations.

<img src="assets/flippy.png" alt="flippy" width="300"/>

[[_TOC_]]

## Gallery

<img src="https://surfdrive.surf.nl/files/index.php/s/6HCtX7B4NwE6w48/download" alt="rbc" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/1O3oJBNMT0vLxIv/download" alt="bubble_collision" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/mua39VAvdxmobD4/download" alt="cell division" width="300">

## Current Version

current release [![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases) is experimental. This means that the API may change significantly. 

## Support
Flippy is still in active development and the documentation is almost non-existent, but I try my best to reply to e-mails and provide support for scientists who want to use flippy.
### for questions about general usage
please use the support email [![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org).
### for bugfixes 
please create an [issue](https://gitlab.tudelft.nl/idema-group/flippy/-/issues).
### for feature requests
again the [issues](https://gitlab.tudelft.nl/idema-group/flippy/-/issues) page can be used but be aware that new features will be slow to come.
### want to join the team?
write me an email!


## Release Version Nomenclature

This repository follows [Semantic Versioning](https://semver.org/) guidelines.

Given a version number *MAJOR*.*MINOR*.*PATCH*, the flippy project will increment the:

- *MAJOR* version when you make incompatible API changes,
- *MINOR* version when you add functionality in a backwards compatible manner, and
- *PATCH* version when you make backwards compatible bug fixes.

Additional labels for pre-release and build metadata are available as extensions to the *MAJOR*.*MINOR*.*PATCH* format.

## licence 

*flippy* is under MIT License, which means you can do pretty much anything with it. For more information read the `LICENCE` file.

## How to get it

*flippy* is a headers only library, so all you need to is to download the `flippy` subfolder and copy it into your project.

## Documentation
You can find flippy's [user manual](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/User_Manual) and automatically generated code [documentation](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/Documentation) over on the [wiki](https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/docs/mainpage.md).

## Examples of usage

This is a simple example of a fluctuating sphere in under 100 lines. The code is not as elegant as it could be, but it is general enough to be used as a starting pint in a real project. This code assumes that the `flippy` project is in the same folder as the main file.
This code saves a data file in the folder where the binary file will be executed, under the name `test_run.json` which is a full snapshot of the final configuration.
```cpp
//demo/simplest_MC/main.cpp

#include "flippy.hpp"
#include <iostream> // needed for saving the the data file
#include <random> // needed for random displacement generation
#include <ctime> // needed to seed the random number generator
#include <vector>

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easily passed to functions
    double bending_rigidity, K_V, K_A, tension, R, V_t, A_t, dx, l_max_square;
    int n_triang, t_max;
};

template<typename Real, typename Index>
Real surface_energy_tension_volume_ensemble(fp::Node<Real, Index> const& node, fp::Triangulation<Real, Index> const& triangulation, SimulationParameters const& prms){
    Real V = triangulation.global_geometry().volume;
    Real A = triangulation.global_geometry().area;
    Real dV = V-prms.V_t;
    Real e_tot = (prms.bending_rigidity/2.)*node.scaled_curvature_energy + prms.K_V*dV*dV/prms.V_t + prms.tension*A;
    return e_tot;
}


template<typename Real, typename Index>
Real surface_energy_area_volume_ensemble(fp::Node<Real, Index> const& node, fp::Triangulation<Real, Index> const& triangulation, SimulationParameters const& prms){
    Real V = triangulation.global_geometry().volume;
    Real A = triangulation.global_geometry().area;
    Real dV = V-prms.V_t;
    Real dA = A-prms.A_t;
    Real e_tot = (prms.bending_rigidity/2.)*node.scaled_curvature_energy + prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return e_tot;
}

int main(){
    fp::print("starting"); // write the string "starting" to the standard out. fp is flippy's namespace and print is a built-in print function
    fp::Timer timer; //setting up a timer that will print to console how long the simulation took
    float R = 10; // initial radius of the triangulated sphere
    int n_triang = 10; // triangulation iteration number of nodes N_node=12+30*n_triang+20*n_triang(n_triang-1)/2
    double l_min = 2.*R*sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.)); // estimate of a typical bond length in the initial triangulation. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 2*l_min;
    const SimulationParameters prms{
        .bending_rigidity = 20 /*kBT_*/, .K_V = 1000 /*kBT_*/, .K_A=1000, .tension = 0.01/(R*R) /*kBT_/R^2*/, // physical parameters that enter the energy
        .R=R /*a.u.*/, .V_t=0.6*sphere_vol(R),
        .A_t=sphere_area(R),
        .dx=l_max/100., // side length of a voxel from which the displacement of the node is drawn
        .l_max_square=l_max*l_max, // max allowed bond length. If bond length of 
        .n_triang=n_triang, 
        .t_max=1000 // max number of iteration steps (depending on the strength of your cpu, this should take anywhere from a couple of seconds to a couple of minutes
    };

    std::random_device random_number_generator_seed;
    std::mt19937 rng(random_number_generator_seed()); // create a random number generator and seed it with current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, int> tr(prms.n_triang, prms.R, l_max);
    fp::MonteCarloUpdater<double, int, SimulationParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION> mc_updater(tr, prms, surface_energy_area_volume_ensemble<double, int>, rng, l_min, l_max);

    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector
    std::uniform_real_distribution<double> dx_distr(-prms.dx, prms.dx); //define a distribution from which the small displacements in x y and z directions will be drawn

    fp::Json data_init = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_init.json", data_init);  // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<int> shuffled_ids;
    shuffled_ids.reserve(tr.size());
    for(auto const& node: tr.nodes()){
        shuffled_ids.push_back(node.id);
    }

    for(int t=0; t<prms.t_max;++t){
        std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
        for (int node_id: shuffled_ids){
            displ={dx_distr(rng), dx_distr(rng), dx_distr(rng)};
            mc_updater.move_MC_updater(tr[node_id], displ);
        }
        for (int node_id: shuffled_ids){
            mc_updater.flip_MC_updater(tr[node_id]);
        }
    }

    // MonteCarloUpdater counts the number of accepted and rejected moves, and it distinguishes between if a rejection
    // occurred because of the energy or the bond length constraint. We can use this to print a simple statistics here.
    // This will help us decide if our displacement size is too large for example.
    fp::print("percentage of failed moves: ",(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((double)mc_updater.move_attempt_count()));
    fp::print("percentage of failed flips: ",(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((double)mc_updater.flip_attempt_count()));

    fp::Json data_final = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_final.json", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable
    timer.stop(); // strictly speaking this is not necessary the timer would stop and print the time automatically when it gets deleted
    return 0;
}

```
This code can be found in the subfolder `demo/simplest_MC` together with a python script that visualizes the data, and a simple `CMake` file that can be used to build the code.
