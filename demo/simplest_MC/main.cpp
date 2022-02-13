#include <random> // needed for random displacement generation
#include <vector>
#include "flippy.hpp"

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct SimulationParameters{ // a data structure that can hold all simulation parameter, s.t. they can be easily passed to functions
    double bending_rigidity, K_V, K_A, R, V_t, A_t, linear_displ, l_max_square;
    int n_triang, max_mc_steps, surface_updates_per_mc_step;
}tension;

// This is the energy function that is used by flippy's built in updater to decide if a move was energetically favourable or not
double surface_energy_area_volume_ensemble([[maybe_unused]]fp::Node<double, int> const& node, fp::Triangulation<double, int> const& triangulation, SimulationParameters const& prms){
    double V = triangulation.global_geometry().volume;
    double A = triangulation.global_geometry().area;
    double dV = V-prms.V_t;
    double dA = A-prms.A_t;
    double e_tot = (prms.bending_rigidity/2.)*triangulation.global_geometry().dA_K2 + prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return e_tot;
}

int main(){
    fp::print("starting"); // write the string "starting" to the standard out. fp is flippy's namespace and print is a built-in print function
    fp::Timer timer; //setting up a timer that will print to console how long the simulation took
    int n_triang = 14; // triangulation iteration number of nodes N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_triang
    double l_min = 2;
    double R = l_min/(2*sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.))); // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length are close to minimal. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 2.5*l_min; // if you make l_max closer to l_min bond_flip acceptance rate will go down
    double V0 = sphere_vol(R);
    SimulationParameters prms{
        .bending_rigidity = 10 /*kBT_*/, .K_V = 100 /*kBT_*/, .K_A=1000, .R=R /*a.u.*/,
        .V_t=V0, .A_t=sphere_area(R),
        .linear_displ=l_min/8., // side length of a voxel from which the displacement of the node is drawn
        .n_triang=n_triang,
        .max_mc_steps=300, // max number of iteration steps (depending on the strength of your cpu, this should take anywhere from a couple of seconds to a couple of minutes
        .surface_updates_per_mc_step=100
    };

    std::random_device random_number_generator_seed;
    std::mt19937 rng(random_number_generator_seed()); // create a random number generator and seed it with current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, int> tr(prms.n_triang, prms.R, 2*l_min);
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

    for(int t=1; t<prms.max_mc_steps+1; ++t){
        if(t%10==0){
            fp::print("V/V_t=",tr.global_geometry().volume/prms.V_t, "A/A_t=", tr.global_geometry().area/prms.A_t);
            tr.make_global_geometry();
            tr.make_verlet_list();
        }
        progress = t/((double)prms.max_mc_steps);
        if(t<prms.max_mc_steps/2){
            //for the first half of the simulation we gradually decrease the target volume. If we do not do this the volume term will dominate the energy, and we will get weird shapes!
            // since the MonteCarloUpdater takes a reference to the SimulationParameters struct, we can just change the parameters and flippy will have access to the updated version.
            prms.V_t = V0*(1.-0.4*(2*progress));
        }

        else if(t<(2*prms.max_mc_steps)/3){
            mc_updater.reset_kBT(1. - (progress-0.5)); //here we cool down the temperature to force the monte carlo simulation to go to the minimum configuration faster
        }
        for(int update_step=0; update_step<prms.surface_updates_per_mc_step;++update_step) {
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
    fp::print("percentage of failed moves: ",(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((double)mc_updater.move_attempt_count()));
    fp::print("percentage of failed flips: ",(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((double)mc_updater.flip_attempt_count()));

    fp::Json data_final = {{"nodes", tr.make_egg_data()}};
    fp::json_dump("test_run_final", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable
    timer.stop(); // strictly speaking this is not necessary the timer would stop and print the time automatically when it gets deleted
    return 0;
}