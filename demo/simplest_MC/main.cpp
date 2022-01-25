#include "flippy/flippy.hpp"
#include "flippy/Triangulator.hpp"
#include <iostream> // needed for saving the the data file
#include <random> // needed for random displacement generation
#include <ctime> // needed to seed the random number generator

float sphere_vol(float R){return 4.f/3.f * M_PI *R*R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easilly passed to functions
    float bending_rigidity, K_V, tension, R, V, dx, l_max_square;
    int n_triang, t_max;
};

template<typename Real, typename Index>
Real surface_energy_tension_volume_ensemble(fp::Node<Real, Index> const& node, fp::Triangulation<Real, Index> const& triangulation, SimulationParameters const& prms){
    Real V = triangulation.global_geometry().volume;
    Real A = triangulation.global_geometry().area;
    Real dV = V-prms.V;
    Real e_tot = (prms.bending_rigidity/2.)*node.scaled_curvature_energy + prms.K_V*dV*dV/V + prms.tension*A;
    return e_tot;
}

template<typename Real>
bool move_needs_undoing(Real e_old, Real e_new, std::mt19937& rng){
    Real e_diff = e_old - e_new;
    return (e_diff<0)&&(std::uniform_real_distribution<Real>(0,1)(rng)>std::exp(e_diff));
} 

template<typename Real, typename Index>
bool next_neighbour_distances_are_smaller_than_max_length(fp::Node<Real, Index> const& node, fp::vec3<Real> const & displacement, Real l_max_square){
    for (auto const& nn_dist: node.nn_distances) // iterate through the next neighbour distances of a node which is a collection of 3d vectors, i.e. nn_dist is of type fp::vec3
    {
        if((nn_dist-displacement).norm_square() > l_max_square){return false;}
    }
    return true;
}

template<typename Real, typename Index>
void move_MC_updater(fp::Node<Real, Index>const& node, fp::vec3<Real> const& displacement, fp::Triangulation<Real, Index> & triangulation, SimulationParameters const& prms, std::mt19937& rng, 
        std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index> const&, SimulationParameters const&)> energy_function){
    if (next_neighbour_distances_are_smaller_than_max_length(node, displacement, prms.l_max_square)){
        Real e_old = energy_function(node, triangulation, prms);
        triangulation.move_node(node.id, displacement);
        Real e_new = energy_function(node, triangulation, prms);
        if(move_needs_undoing<Real>(e_old, e_new, rng)){triangulation.move_node(node.id,-displacement);}
    }
}

template<typename Real, typename Index>
void flip_MC_updater(fp::Node<Real, Index>const& node, fp::Triangulation<Real, Index> & triangulation, SimulationParameters const& prms, std::mt19937& rng, 
        std::function<Real(fp::Node<Real, Index> const&, fp::Triangulation<Real, Index> const&, SimulationParameters const&)> energy_function){
    Real e_old = energy_function(node, triangulation, prms);
    int num_nns = node.nn_ids.size();
    int nn_id = node.nn_ids[rand()%num_nns];
    auto bfd = triangulation.flip_bond(node.id, nn_id, 0., prms.l_max_square);
    if(bfd.flipped){
        Real e_new = energy_function(node, triangulation, prms);
        if(move_needs_undoing(e_old, e_new, rng)){triangulation.unflip_bond(node.id, nn_id, bfd);}
    }
}

int main(){
    fp::print("starting"); // write the string "starting" to the standard out. fp is flippy's namespace and print is a built in print function
    fp::Timer timer; //setting up a timer that will print to console how long the simulation took
    float R = 10; // initial radius of the triangulated sphere
    int n_triang = 14; // triangulation iteration
    float l_bond_init_est = R/(n_triang+1.f); // estimate of a typical bond length in the initial triangulation. Tis varies a bit around R/(nIt+1) due to implementation details
    float l_max = 1.4f*l_bond_init_est;
    const SimulationParameters prms{
        .bending_rigidity = 20 /*kBT*/, .K_V = 100 /*kBT*/, .tension = 0.01f/(R*R) /*kBT/R^2*/, // physical parameters that enter the energy
        .R=R /*a.u.*/, .V=sphere_vol(R), 
        .dx=l_max/5.f, // side leng of a voxel from which the displacement of the node is drawn
        .l_max_square=l_max*l_max, // max allowed bond length. If bond length of 
        .n_triang=n_triang, 
        .t_max=3000 // max number of iteration steps
    };
    
    std::mt19937 rng(time(0)); // create a random number generator and seed it with current time 
    fp::Triangulation<float, int> tr(prms.n_triang, prms.R, 1.f);
    fp::vec3<float> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector 
    std::uniform_real_distribution<float> dx_distr(-prms.dx, prms.dx); //define a distribution from which the small displacements in x y and z directions will be drawn

    for(int t=0; t<prms.t_max;++t){
        for (auto const& node: tr.nodes()){
            displ={dx_distr(rng), dx_distr(rng), dx_distr(rng)};
            move_MC_updater<float, int>(node, displ, tr, prms, rng, surface_energy_tension_volume_ensemble<float,int>);
            flip_MC_updater<float, int>(node, tr, prms, rng, surface_energy_tension_volume_ensemble<float,int>);
        }
    }

    fp::Json data = {{"nodes", tr.egg_data()}};
    std::ofstream file("test_run.json");
    file<<data;
    timer.stop(); // strictly speaking this is not necessary the timer would stop and print the time automatically when it gets deleted
    return 0;
}
