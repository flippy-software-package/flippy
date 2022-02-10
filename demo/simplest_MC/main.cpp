#include "flippy.hpp"
#include "Triangulator.hpp"
#include <iostream> // needed for saving the the data file
#include <random> // needed for random displacement generation
#include <ctime> // needed to seed the random number generator
#include <vector>

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easilly passed to functions
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
    fp::print("starting"); // write the string "starting" to the standard out. fp is flippy's namespace and print is a built in print function
    fp::Timer timer; //setting up a timer that will print to console how long the simulation took
    float R = 10; // initial radius of the triangulated sphere
    int n_triang = 12; // triangulation iteration number of nodes N_node=12+30*n_triang+20*n_triang(n_triang-1)/2
    double l_min = 2.*R*sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.)); // estimate of a typical bond length in the initial triangulation. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 1.4*l_min;
    const SimulationParameters prms{
        .bending_rigidity = 20 /*kBT*/, .K_V = 100 /*kBT*/, .K_A=100, .tension = 0.01/(R*R) /*kBT/R^2*/, // physical parameters that enter the energy
        .R=R /*a.u.*/, .V_t=0.6*sphere_vol(R),
        .A_t=sphere_area(R),
        .dx=l_max/10., // side leng of a voxel from which the displacement of the node is drawn
        .l_max_square=l_max*l_max, // max allowed bond length. If bond length of 
        .n_triang=n_triang, 
        .t_max=3000 // max number of iteration steps
    };
    
    std::mt19937 rng(time(0)); // create a random number generator and seed it with current time 
    fp::Triangulation<double, int> tr(prms.n_triang, prms.R, l_max);
    fp::MonteCarloUpdater<double, int, SimulationParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION> mc_updater(tr, prms, surface_energy_area_volume_ensemble<double, int>, rng, l_min, l_max);
    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector
    std::uniform_real_distribution<double> dx_distr(-prms.dx, prms.dx); //define a distribution from which the small displacements in x y and z directions will be drawn

    fp::Json data_init = {{"nodes", tr.make_egg_data()}};
    std::ofstream file_init("test_run_init.json");
    file_init<<data_init;

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

    fp::Json data = {{"nodes", tr.make_egg_data()}};
    std::ofstream file("test_run.json");
    file<<data;
    timer.stop(); // strictly speaking this is not necessary the timer would stop and print the time automatically when it gets deleted
    return 0;
}
