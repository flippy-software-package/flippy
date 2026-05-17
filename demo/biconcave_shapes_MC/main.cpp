#include "flippy.hpp"
#include <iostream> // needed for std::cout
#include <random>   // needed for random displacement generation and node index shuffling
#include <vector>   // need for std::vector

namespace {
using fp::operator""_r;

fp::Real sphere_vol(fp::Real R) { return (4_r / 3_r) * fp::PI * R * R * R; }
fp::Real sphere_area(fp::Real R) { return 4_r * fp::PI * R * R; }

struct EnergyParameters {
    fp::Real kappa, K_V, K_A, V_t, A_t;
};

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable
// or not
fp::Real surface_energy(const fp::Node               &node,
                        const fp::Triangulation      &trg,
                        const EnergyParameters       &p,
                        const std::vector<fp::Index> &changed_neighbourhood) {
    fp::Real V         = trg.global_geometry().volume;
    fp::Real A         = trg.global_geometry().area;
    fp::Real dV        = V - p.V_t;
    fp::Real dA        = A - p.A_t;
    fp::Real E_bending = fp::node_unit_bending_energy(node);
    for (fp::Index changed_node_id : changed_neighbourhood) {
        E_bending += fp::node_unit_bending_energy(trg[changed_node_id]);
    }

    fp::Real energy = (p.kappa * E_bending) + (p.K_V * dV * dV / p.V_t) + (p.K_A * dA * dA / p.A_t);
    return energy;
}

} // namespace

int main() {
    // triangulation iteration number of nodes N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
    //
    fp::Index n_triang = 7;
    fp::Real  l_min    = 2_r;

    // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial
    // bond length is close to minimal. This formula is derived from the equidistant sub-triangulation of an
    // icosahedron, where geodesic distances are used as a distance measure.
    fp::Real R = fp::min_radius_with_non_overlapping_beads(l_min, n_triang);

    fp::Real l_max    = 1.8 * l_min; // if you make l_max closer to l_min bond_flip acceptance rate will go down
    fp::Real r_Verlet = 2 * l_max;
    fp::Real red_vol  = 0.6;

    EnergyParameters prms{.kappa = 10_r /*kBT*/,
                          .K_V   = 100_r /*kBT/area*/,
                          .K_A   = 1000_r /*kBT/volume*/,
                          .V_t   = red_vol * sphere_vol(R),
                          .A_t   = sphere_area(R)};
    fp::Real linear_displ = l_min / 18_r; // side length of a voxel from which the displacement of the node is drawn
    int max_mc_steps = 2e5; // max number of iteration steps (depending on the strength of your CPU, this should take
                            // anywhere from a couple of seconds to a couple of minutes

    std::random_device random_number_generator_seed;
    std::mt19937       rng(
        random_number_generator_seed()); // create a random number generator and seed it with the current time

    // All the flippy magic is happening on the following two lines
    auto guv = fp::Triangulation::make_spherical_triangulation(n_triang, R, r_Verlet);
    fp::MonteCarloUpdater<EnergyParameters, std::mt19937> mc_updater(guv, prms, surface_energy, rng, l_min, l_max);

    fp::vec3<fp::Real> displ{}; // declaring a 3d vector (using flippy's built in vec3 type) for later use as a random
    // direction vector
    std::uniform_real_distribution<fp::Real> displ_distr(-linear_displ, linear_displ); // define a distribution from
    // which the small displacements in x y and z directions will be drawn

    guv.scale_node_coordinates(1, 1, 0.85); // squish the sphere in the
    // z-direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly.

    fp::Json data_init = guv.make_egg_data();
    fp::json_dump("test_run_init",
                  data_init); // ATTENTION!!! this file will be saved in the same folder as the executable

    std::vector<fp::Index> shuffled_ids;
    shuffled_ids.reserve(guv.size());
    for (const auto &node : guv.nodes()) {
        shuffled_ids.push_back(node.id);
    } // create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly
      // through the nodes

    for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
        for (fp::Index node_id : shuffled_ids) { // we first loop through all the beads and move them
            displ = {.x = displ_distr(rng), .y = displ_distr(rng), .z = displ_distr(rng)};
            mc_updater.move_MC_updater(guv[node_id], displ); // guv[node_id] returns the node which has id=node_id
        }
        std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng); // then we shuffle the bead_ids
        for (fp::Index node_id : shuffled_ids) { // then we loop through all of them again and try to flip their bonds
            mc_updater.flip_MC_updater(guv[node_id]);
        }
    }

    // MonteCarloUpdater counts the number of accepted and rejected moves, distinguishing whether a rejection occurred
    // because of the energy or the bond length constraint. We can use this to print simple statistics here. For
    // example, this will help us decide if our displacement size is too large.
    std::cout << "percentage of failed moves: "
              << static_cast<fp::Real>(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count()) /
                     static_cast<fp::Real>(mc_updater.move_attempt_count())
              << '\n';
    std::cout << "percentage of failed flips: "
              << static_cast<fp::Real>(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count()) /
                     static_cast<fp::Real>(mc_updater.flip_attempt_count())
              << '\n';

    fp::Json data_final = guv.make_egg_data();
    fp::json_dump("test_run_final",
                  data_final); // ATTENTION!!! this file will be saved in the same folder as the executable

    return 0;
}
