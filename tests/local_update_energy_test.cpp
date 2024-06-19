#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "flippy.hpp"

using REAL = float;
using INDEX = unsigned int;
using Catch::Matchers::WithinAbs;
[[maybe_unused]] static auto kinda_close(auto num){
    return (WithinAbs(static_cast<REAL>(num),
                      static_cast<REAL>(0.02)));
}

struct EnergyParameters{REAL kappa, K_V, K_A, V_t, A_t;};

// This is the energy function that is used by flippy's built-in updater to decide if a move was energetically favorable or not
static REAL surface_energy([[maybe_unused]]fp::Node<REAL, INDEX> const& node,
                    fp::Triangulation<REAL, INDEX> const& trg,
                    EnergyParameters const& prms, std::optional<std::array<INDEX, 3>> flip_partners={}){

    REAL eb_local = node.unit_bending_energy;

    for (auto const& nn_id: node.nn_ids) {
        eb_local += trg[nn_id].unit_bending_energy;
    }

    if(flip_partners.has_value()){
        std::vector<INDEX> updated_neighborhood;
        updated_neighborhood.reserve(30);
        updated_neighborhood = node.nn_ids;
        updated_neighborhood.push_back(node.id);
        for (auto const& flip_partner_id: flip_partners.value()) {
//            if(!fp::is_member(unique_neighbours, flip_partner_id)){ unique_neighbours.push_back(flip_partner_id);}
            for(auto const& nn_id: trg[flip_partner_id].nn_ids)
                if(!fp::is_member(updated_neighborhood, nn_id)){
                    updated_neighborhood.push_back(nn_id);
                }
        }

        for(auto const& node_id: updated_neighborhood){
            eb_local += trg[node_id].unit_bending_energy;
        }

    }

    REAL V = trg.global_geometry().volume;
    REAL A = trg.global_geometry().area;

    REAL dV = V-prms.V_t;
    REAL dA = A-prms.A_t;

    REAL energy = prms.kappa*eb_local +
                  prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
//    REAL energy = prms.kappa*trg.global_geometry().unit_bending_energy +
//                  prms.K_V*dV*dV/prms.V_t + prms.K_A*dA*dA/prms.A_t;
    return energy;
}

TEST_CASE("local update test"){
    using Triangulation = fp::Triangulation<REAL, INDEX, fp::SPHERICAL_TRIANGULATION>;
    using MCU = fp::MonteCarloUpdater<REAL, INDEX, EnergyParameters, std::mt19937, fp::SPHERICAL_TRIANGULATION>;
    // triangulation iteration number of nodes
    // N_node=12+30*n+20*n*(n-1)/2 where n is the same as n_trng
    uint n_triang = 8;
    REAL l_min = 2;
    REAL kappa = 10; /*kBT*/
    REAL K_A =  1000;/*kBT/volume*/
    REAL K_V = 1000; /*kBT/area*/
    REAL red_vol = 0.6f;
    int max_mc_steps = 1e5;
    std::string save_dir = "../../demo_out/local_update_test/";


    // estimate of a typical bond length in the initial triangulation and then create a sphere such that the initial bond length is close to minimal. This formula is derived from the equidistant sub-triangulation of an icosahedron, where geodesic distances are used as a distance measure.
    REAL R = fp::min_radius_With_non_overlapping_beads(l_min, n_triang);
    REAL l_max = 1.7f * l_min; // if you make l_max closer to l_min
    // bond_flip acceptance rate will go down
    REAL r_Verlet = 2.f * l_max;

    REAL Vt = red_vol * fp::sphere_vol(R);
    REAL At = fp::sphere_area(R);
    EnergyParameters prms{.kappa= kappa, .K_V=K_V, .K_A=K_A, .V_t=Vt, .A_t=At};
    // side length of a voxel from which the displacement of the node is drawn
    REAL linear_displ = l_min / 12.f;
    // max number of iteration steps (depending on the strength of your CPU, this should take anywhere from a couple of seconds to a couple of minutes

    std::random_device random_number_generator_seed;
    // create a random number generator and seed it with the current time
    std::mt19937 rng(random_number_generator_seed());

    auto guv = Triangulation(n_triang, R, r_Verlet);
    auto mc_updater = MCU(guv, prms, surface_energy, rng, l_min, l_max);

    fp::vec3<REAL> displ{};
    std::uniform_real_distribution<REAL> displ_distr(-linear_displ, linear_displ);

    // squish the sphere in the z-direction to break the initial symmetry. This speeds up the convergence to a biconcave shape greatly.
//    guv.scale_node_coordinates(0.90f, 0.90f, 1.f);
    guv.scale_node_coordinates(1.f, 1.f, 0.8f);

    fp::Json data_init = guv.make_egg_data();
    fp::json_dump(save_dir+"test_run_init", data_init);

    std::vector<INDEX> shuffled_ids;
    shuffled_ids.reserve(guv.size());

    //create a vector that contains all node ids. We can shuffle this vector in each MC step to iterate randomly through the nodes
    for (auto const &node: guv.nodes()) {
        shuffled_ids.push_back(node.id);
    }
    auto measured_porob = [](MCU& mcu){
        double p_rej = static_cast<double>(mcu.move_back_count() + mcu.bond_length_move_rejection_count()) / static_cast<double>(mcu.move_attempt_count());
        return static_cast<REAL>(1.-p_rej);
    };

    SECTION("always reject: ") {
        REAL probability_target = 0.4f;
        REAL adaptation_magnitude = 0.1f;
        std::optional<REAL> e_diff = 0;

        REAL V0 = guv.global_geometry().volume;
        REAL A0 = guv.global_geometry().area;
        auto displ_updater = fp::DynamicDisplacementUpdater<REAL, INDEX>(linear_displ, adaptation_magnitude, probability_target);

        for (int mc_step = 0; mc_step < max_mc_steps; ++mc_step) {
            auto t = static_cast<REAL>(mc_step);
            auto t_max = static_cast<REAL>(max_mc_steps);
            prms.V_t = fp::linear_adaptation(t, 0.f, 0.3f*t_max, V0, Vt);
            prms.A_t = fp::linear_adaptation(t, 0.f, 0.3f*t_max, A0, At);
            REAL kT = fp::linear_adaptation(t, 0.7f*t_max, 0.9f*t_max, 1.f, 0.1f);
            mc_updater.reset_kBT(kT);
            displ_distr = displ_updater.new_displ_distr();
            for (INDEX node_id: shuffled_ids) {
                displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
                e_diff = mc_updater.move_MC_updater(guv[node_id], displ);
                displ_updater.probability_aggregator(e_diff, mc_updater.kBT());
            }

            REAL pt = fp::linear_adaptation(t, 0.7f*t_max, 0.9f*t_max, probability_target, 0.2f);
            displ_updater.reset_target_acceptance_probability(pt);
            std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng);
            for (auto node_id: shuffled_ids) { mc_updater.flip_MC_updater(guv[node_id]); }
        }
        CHECK_THAT(measured_porob(mc_updater), kinda_close(probability_target));
        std::cout<<"rejected moves: "<<static_cast<long double>(mc_updater.move_back_count() + mc_updater.bond_length_move_rejection_count()) / static_cast<long double>(mc_updater.move_attempt_count() ) <<'\n';
        std::cout<<"rejected flips: "<<static_cast<long double>(mc_updater.flip_back_count() + mc_updater.bond_length_flip_rejection_count()) / static_cast<long double>(mc_updater.flip_attempt_count()) <<'\n';
    }

    fp::Json data_final = guv.make_egg_data();
    fp::json_dump(save_dir+"test_run_final", data_final);


}
