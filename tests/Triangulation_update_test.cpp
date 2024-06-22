#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <array>
#include <iostream>

#include "flippy.hpp"
#include "catch2/generators/catch_generators_all.hpp"

using namespace fp;

//template <floating_point_number Real, indexing_number Index>
//void rescale_triangulation(Real R, Triangulation<Real,Index, SPHERICAL_TRIANGULATION>& tr)
//{
//    tr.R_initial=R;
//    tr.scale_all_nodes_to_R_init();
//    tr.orient_surface_of_a_sphere();
//    tr.initiate_distance_vectors(); //Todo if this is done before orient surface we can save time
//    tr.make_global_geometry();
//    tr.make_verlet_list();
//}

bool euler_number_is_2(Triangulation<SPHERICAL_TRIANGULATION>const& triangulation)
{
    /**
     * General euler formula for flat polyhedra: `V - E + F = 2`
     * Special formula for triangulations is derived through `F = (2/3)E`:
     *
     * `E = 3V - 6`
     *
     */
    Index V = triangulation.size();
    Index E = 0;

    for(auto const& node : triangulation.nodes()){ E += node.nn_ids.size(); }

    E = E/2;

    return E == 3*V - 6;

}

SCENARIO("Random bond flips do not destroy the topology"){
    GIVEN("A a spherical triangulation"){
        constexpr Index subriang_level = 10;
        constexpr Real verlet_radius = 0.f;
        constexpr Real r_init = 13.f;
        constexpr Index test_repetitions = 10;
        Triangulation<SPHERICAL_TRIANGULATION> tr(subriang_level, r_init, verlet_radius);
        auto rg = Catch::Generators::RandomIntegerGenerator<Index>(0, tr.size(), Catch::Generators::Detail::getSeed());

        THEN("Euler Number is 2"){
            CHECK(euler_number_is_2(tr));
        }

        WHEN("Bonds are flipped randomly"){
            Index rn = rg.get();
            rg.get();

            for(Index _=0; _<test_repetitions;++_) {
                for (auto const &node: tr.nodes()) {
                    auto n_neighbors = static_cast<Index>(node.nn_ids.size());
                    Index rand_id = rn % n_neighbors;
                    tr.flip_bond(node.id, node.nn_ids[rand_id], 0, r_init*r_init);
                }
            }

            THEN("The the euler characteristic of the triangulation does not change"){
                CHECK(euler_number_is_2(tr));
            }

        }
    }
}

