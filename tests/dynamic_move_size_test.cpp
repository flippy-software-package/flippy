#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "flippy.hpp"


[[maybe_unused]] static auto kinda_abs_close(auto num){
    using Catch::Matchers::WithinAbs;
    if(num*num < 10e-10){}
    return (WithinAbs(static_cast<fp::Real>(num),
                      static_cast<fp::Real>(0.01)));
}

[[maybe_unused]] static auto kinda_rel_close(auto num){
    using Catch::Matchers::WithinRel;
    if(num*num < 10e-10){}
    return (WithinRel(static_cast<fp::Real>(num),
                      static_cast<fp::Real>(0.05)));
}

struct EnergyParameters{fp::Real kappa, K_V, K_A, V_t, A_t;};
using Triangulation = fp::Triangulation;
using MCU = fp::MonteCarloUpdater<EnergyParameters, std::mt19937>;

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable or not
fp::Real surface_energy([[maybe_unused]]fp::Node const& node,
                    fp::Triangulation const& trg,
                    EnergyParameters const& prms,
                    std::vector<fp::Index>const& changed_neighborhood){
    fp::Real V = trg.global_geometry().volume;
    fp::Real A = trg.global_geometry().area;
    fp::Real dV = V-prms.V_t;
    fp::Real dA = A-prms.A_t;

    fp::Real eb = fp::changed_neighborhood_bending_energy(node, trg.nodes(), changed_neighborhood);

    fp::Real energy = prms.kappa*eb +
                  prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return energy;
}

TEST_CASE("spherical triangulation, with initially broken symmetry (stretched to an ellipse)"){
        // triangulation iteration number of nodes
        // N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
        uint n_triang = 4;
        fp::Real l_min = 2;
        fp::Real kappa = 10; /*kBT*/
        fp::Real K_A =  1000;/*kBT/volume*/
        fp::Real K_V = 100; /*kBT/area*/
        fp::Real red_vol = 0.64f;
        int max_mc_steps = 7e3;
        std::string save_dir = "./";


        // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length is close to minimal. This formula is derived from the equidistant sub-triangulation of an icosahedron, where geodesic distances are used as a distance measure.
        fp::Real R = fp::min_radius_with_non_overlapping_beads(l_min, n_triang);
        fp::Real l_max = 2.f * l_min; // if you make l_max closer to l_min
        // bond_flip acceptance rate will go down
        fp::Real r_Verlet = 2.f * l_max;

        EnergyParameters prms{.kappa= kappa, .K_V=K_V, .K_A=K_A, .V_t=red_vol * fp::sphere_vol(R), .A_t=fp::sphere_area(R)};
        // side length of a voxel from which the displacement of the node is drawn
        fp::Real linear_displ = l_min / 8.f;
        // max number of iteration steps (depending on the strength of your CPU, this should take anywhere from a couple of seconds to a couple of minutes

        std::random_device random_number_generator_seed;
        // create a random number generator and seed it with the current time
        std::mt19937 rng(random_number_generator_seed());

        auto guv = Triangulation(n_triang, R, r_Verlet);
        auto mc_updater = MCU(guv, prms, surface_energy, rng, l_min, l_max);

        fp::vec3<fp::Real> displ{};
        std::uniform_real_distribution<fp::Real> displ_distr(-linear_displ, linear_displ);

        // squish the sphere in the z-direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly.
        guv.scale_node_coordinates(1.f, 1.f, 0.8f);

        fp::Json data_init = guv.make_egg_data();
        fp::json_dump(save_dir+"test_run_init", data_init);

        std::vector<fp::Index> shuffled_ids;
        shuffled_ids.reserve(guv.size());

        //create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly through the nodes
        for (auto const &node: guv.nodes()) {
            shuffled_ids.push_back(node.id);
        }
    auto measured_porob = [](MCU& mcu){
        double p_rej = static_cast<double>(mcu.move_back_count() + mcu.bond_length_move_rejection_count()) / static_cast<double>(mcu.move_attempt_count());
        return static_cast<fp::Real>(1.-p_rej);
        };

    SECTION("always accept: ") {
        fp::Real probability_target = 1.f;
        std::optional<fp::Real> e_diff = 0;

        auto displ_updater = fp::DynamicDisplacementUpdater(linear_displ, probability_target);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            displ_distr = displ_updater.new_displ_distr();
            for (fp::Index node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (auto node_id: shuffled_ids) { mc_updater.flip_MC_updater(guv[node_id]); }
        }
        CHECK_THAT(measured_porob(mc_updater), kinda_rel_close(probability_target));
    }


    SECTION("always reject: ") {
        fp::Real probability_target = 0.f;
        std::optional<fp::Real> e_diff = 0;

        auto displ_updater = fp::DynamicDisplacementUpdater(linear_displ, probability_target);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            displ_distr = displ_updater.new_displ_distr();
            for (fp::Index node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (auto node_id: shuffled_ids) { mc_updater.flip_MC_updater(guv[node_id]); }
        }
        CHECK_THAT(measured_porob(mc_updater), kinda_abs_close(probability_target));
    }

    SECTION("fifty fifty: ") {
        fp::Real probability_target = 0.5f;
        std::optional<fp::Real> e_diff = 0;

        auto displ_updater = fp::DynamicDisplacementUpdater(linear_displ, probability_target);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            displ_distr = displ_updater.new_displ_distr();
            for (fp::Index node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (auto node_id: shuffled_ids) { mc_updater.flip_MC_updater(guv[node_id]); }
        }
        CHECK_THAT(measured_porob(mc_updater), kinda_rel_close(probability_target));
    }


    SECTION("sixty percent accept: ") {
        fp::Real probability_target = 0.6f;
        std::optional<fp::Real> e_diff = 0;

        auto displ_updater = fp::DynamicDisplacementUpdater(linear_displ, probability_target);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            displ_distr = displ_updater.new_displ_distr();
            for (fp::Index node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (auto node_id: shuffled_ids) { mc_updater.flip_MC_updater(guv[node_id]); }
        }
        CHECK_THAT(measured_porob(mc_updater), kinda_rel_close(probability_target));
    }

}
