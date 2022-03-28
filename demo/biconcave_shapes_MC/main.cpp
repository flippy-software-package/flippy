#include <random> // needed for random displacement generation
#include <vector> // need for std::vector
#include <iostream> // needed for std::cout
#include "flippy.hpp"

double sphere_vol(double R){return 4./3. * M_PI *R*R*R;}
double sphere_area(double R){return 4. * M_PI *R*R;}

struct EnergyParameters{double kappa, K_V, K_A, V_t, A_t;};

// This is the energy function that is used by flippy's built in updater to decide if a move was energetically favourable or not
double surface_energy([[maybe_unused]]fp::Node<double, int> const& node,
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

int main(){
    int n_triang = 7; // triangulation iteration number of nodes N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
    double l_min = 2;
    double R = l_min/(2*sin(asin(1./(2*sin(2.*M_PI/5.)))/(n_triang+1.))); // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length are close to minimal. This formula is derived from equidistant subtriangulation of an icosahedron, where geodesic distances are used as a distance measure.
    double l_max = 2.*l_min; // if you make l_max closer to l_min bond_flip acceptance rate will go down
    double r_Verlet = 2*l_max;
    EnergyParameters prms{.kappa=10 /*kBT*/,
                          .K_V=100 /*kBT/area*/, .K_A=1000 /*kBT/volume*/,
                          .V_t=0.6*sphere_vol(R), .A_t=sphere_area(R)};
    double linear_displ=l_min/8.; // side length of a voxel from which the displacement of the node is drawn
    int max_mc_steps=1e5; // max number of iteration steps (depending on the strength of your cpu, this should take anywhere from a couple of seconds to a couple of minutes

    std::random_device random_number_generator_seed;
    std::mt19937 rng(random_number_generator_seed()); // create a random number generator and seed it with current time

    // All the flippy magic is happening on the following two lines
    fp::Triangulation<double, int> guv(n_triang, R, r_Verlet);
    fp::MonteCarloUpdater<double, int, EnergyParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION> mc_updater(guv, prms, surface_energy, rng, l_min, l_max);

    fp::vec3<double> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for a later use as a random direction vector
    std::uniform_real_distribution<double> displ_distr(-linear_displ, linear_displ); //define a distribution from which the small displacements in x y and z directions will be drawn

    guv.scale_node_coordinates(1, 1, 0.8); // squish the sphere in z direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly

    fp::Json data_init = guv.make_egg_data();
    fp::json_dump("test_run_init", data_init);  // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<int> shuffled_ids;
    shuffled_ids.reserve(guv.size());
    for(auto const& node: guv.nodes()){ shuffled_ids.push_back(node.id);} //create a vector that contains all node ids. We can shuffle this vector in each MC step, to iterate randomly through the nodes

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

    // MonteCarloUpdater counts the number of accepted and rejected moves, and it distinguishes between if a rejectio occurred because of the energy or the bond length constraint.
    // We can use this to print a simple statistics here. This will help us decide if our displacement size is too large for example.
    std::cout<<"percentage of failed moves: "<<(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count())/((long double)mc_updater.move_attempt_count())<<'\n';
    std::cout<<"percentage of failed flips: "<<(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count())/((long double)mc_updater.flip_attempt_count())<<'\n';

    fp::Json data_final = guv.make_egg_data();
    fp::json_dump("test_run_final", data_final);  // ATTENTION!!! this file will be saved in the same folder as the executable

    return 0;
}